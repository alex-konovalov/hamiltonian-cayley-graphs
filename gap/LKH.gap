# LKH.gap

# This file provides an interface from GAP to LKH (Keld Helsgaun's implementation of the 
# Lin-Kernighan-Helsgaun heuristic). More precisely, the LKH function asks LKH to look
# for a hamiltonian cycle in the (grape) graph X. The program LKH (version 2.0.7) must be 
# installed and be available to GAP. It can be downloaded for academic and noncommercial
# use from
# 	http://www.akira.ruc.dk/~keld/research/LKH/
#
# The function LKH returns the hamiltonian cycle found by LKH, as a list of vertex names 
# in the order that they are visited by the cycle.  It returns fail if LKH did not 
# return a valid hamiltonian cycle that contains all of the required edges.
#
# Inputs are:
# 	X: a GRAPE graph
#	AdditionalEdges: a list of edges (2-element lists) to be added to the graph X.
#	RequiredEdges: a list of edges that are required to be used in the hamiltonian cycle.
# AdditionalEdges and RequiredEdges may be the empty list []. They specify vertices by
# using vertex names. Edges in both lists are assumed to be directed. To specify an 
# undirected edge, also include the same edge with the opposite orientation.

LoadPackage("grape");
Read("IsHamiltonianCycle.gap");

LKHTempDir := DirectoryTemporary();
if LKHTempDir = fail then
	Error("LKH was unable to create a temporary directory for its scratch files\n");
fi;

LKH := function(X, AdditionalEdges, RequiredEdges)

	local LKHtspfile, LKHparfile, LKHtoGAPfile, LKHstdoutfile, LKHstdoutStream,
	LKHtoGAP, AdjMat, edge, StartOfEdge, EndOfEdge, HamCyc, 
	WeightOfNonEdge, WeightOfEdge, WeightOfRequiredEdge,
	i, j;

	WeightOfEdge := 1;
	WeightOfRequiredEdge := 0;
	WeightOfNonEdge := 10;
	
	# a graph with less than 3 vertices cannot have a hamiltonian cycle:
	if OrderGraph(X) <= 2 then return fail; fi; 

	# create files to communicate with the LKH program
    LKHparfile := Filename(LKHTempDir, "GAPtoLKH.par");
    LKHtspfile := Filename(LKHTempDir, "GAPtoLKH.tsp");
    LKHtoGAPfile := Filename(LKHTempDir, "LKHtoGAP.txt");
    LKHstdoutfile := Filename(LKHTempDir, "LKHstdout.txt");
    
    PrintTo(LKHparfile, 
        "PROBLEM_FILE = ", LKHtspfile, "\n",
        "TOUR_FILE = ", LKHtoGAPfile, "\n"
        );

	# if possible, we send an undirected graph to LKH, instead of a directed graph
	if 
		IsSimpleGraph(X)
		and ForAll(AdditionalEdges, edge -> (Reversed(edge) in AdditionalEdges)) 
		and ForAll(RequiredEdges, edge -> (Reversed(edge) in RequiredEdges))
	then
 
		PrintTo(LKHtspfile, 
			"TYPE : HCP\n", # hamiltonian cycle problem (undirected graph)
			"DIMENSION : ", OrderGraph(X), "\n",
			"EDGE_DATA_FORMAT : EDGE_LIST\n",
			"EDGE_DATA_SECTION\n"                # title of section
			);
		
		for edge in UndirectedEdges(X) do
			AppendTo(LKHtspfile, edge[1], " ", edge[2], "\n");
		od;		
		for edge in AdditionalEdges do
			StartOfEdge := Position(VertexNames(X), edge[1]);  # convert start of edge to an integer
			EndOfEdge := Position(VertexNames(X), edge[2]);    # convert end of edge to an integer
			if StartOfEdge < EndOfEdge then
				AppendTo(LKHtspfile, StartOfEdge, " ", EndOfEdge, "\n");
			fi;
		od;
		AppendTo(LKHtspfile, "-1\n");    # -1 marks the end of this section of the data
	
		if not(IsEmpty(RequiredEdges)) then
			AppendTo(LKHtspfile, "FIXED_EDGES_SECTION :\n");
			for edge in RequiredEdges do
				StartOfEdge := Position(VertexNames(X), edge[1]);  # convert start of edge to an integer
				EndOfEdge := Position(VertexNames(X), edge[2]);    # convert end of edge to an integer
				AppendTo(LKHtspfile, StartOfEdge, " ", EndOfEdge, "\n");
			od;    
			AppendTo(LKHtspfile, "-1\n");    # -1 marks the end of this section of the data
		fi;

	else # the graph is directed

		PrintTo(LKHtspfile, 
			"TYPE : ATSP\n", # asymmetric traveling salesman problem (directed graph)
			"DIMENSION : ", OrderGraph(X), "\n",
			"EDGE_WEIGHT_TYPE: EXPLICIT\n",
			"EDGE_WEIGHT_FORMAT: FULL_MATRIX\n",
			"EDGE_WEIGHT_SECTION\n"                # title of section
			);

		# create the adjacency matrix
		# first, initialize to be all non-edges		
		AdjMat := List([1..OrderGraph(X)], i ->  List([1..OrderGraph(X)], j -> WeightOfNonEdge ));
		# put in the edges
		for i in [1..OrderGraph(X)] do
			for j in Adjacency(X, i) do
				AdjMat[i][j] := WeightOfEdge;
			od;
		od;
		# put in the additional edges
		for edge in AdditionalEdges do
			StartOfEdge := Position(VertexNames(X), edge[1]);  # convert start of edge to an integer
			EndOfEdge := Position(VertexNames(X), edge[2]);    # convert end of edge to an integer
			AdjMat[StartOfEdge][EndOfEdge] := WeightOfEdge;
		od;      
		# put in the required edges
		for edge in RequiredEdges do
			StartOfEdge := Position(VertexNames(X), edge[1]);  # convert start of edge to an integer
			EndOfEdge := Position(VertexNames(X), edge[2]);    # convert end of edge to an integer
			AdjMat[StartOfEdge][EndOfEdge] := WeightOfRequiredEdge;
		od;
		# write the adjacency matrix to the tsp file
		for i in [1..OrderGraph(X)] do
			for j in [1..OrderGraph(X)] do
				AppendTo(LKHtspfile, AdjMat[i][j], "\n");
			od;
		od; 
	fi;

	# ask LKH (Keld Helsgaun's Lin-Kernighan-Helsgaun algorithm) to find a hamiltonian cycle
	LKHstdoutStream := OutputTextFile(LKHstdoutfile, false);
	Process(
		LKHTempDir, 
		Filename( DirectoriesSystemPrograms(), "LKH"), 
		InputTextUser(), 
		LKHstdoutStream,
		[LKHparfile]
		);
	CloseStream(LKHstdoutStream);

	# read in the tour that was found by LKH
	LKHtoGAP := InputTextFile(LKHtoGAPfile);
	HamCyc := ReadAll(LKHtoGAP); # read in the raw data (a string)
	CloseStream(LKHtoGAP);
 	
 	# convert the the tour into a GAP list of vertex names
	HamCyc := SplitString(HamCyc, "\n" ); # parse the data into a list of integers
	Remove(HamCyc,1);      # discard the first line (NAME) 
	Remove(HamCyc,1);      # discard the second line (COMMENT - length) 
	Remove(HamCyc,1);      # discard the third line (COMMENT - found by) 
	Remove(HamCyc,1);      # discard the fourth line (TYPE) 
	Remove(HamCyc,1);      # discard the fifth line (DIMENSION) 
	Remove(HamCyc,1);      # discard the sixth line (TOUR_SECTION) 
	Remove(HamCyc);        # discard the last line (EOF) 
	Remove(HamCyc);        # discard the second-to-last line (-1) 
	Add(HamCyc,HamCyc[1]); # place a duplicate of the first vertex at the end (to close the cycle)
	Apply(HamCyc, Int);    # replace each string with the integer it represents
	Apply(HamCyc, x -> VertexName(X, x)); # translate the integers into names of vertices   

	# if HamCyc is a valid hamiltonian cycle, then we return it
	if IsHamiltonianCycle(X, HamCyc, AdditionalEdges, RequiredEdges) then
		return HamCyc;
	else
		return fail;
	fi;
end;

# This is free and unencumbered software released into the public domain.
# 
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
# 
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
# 
# For more information, please refer to <http://unlicense.org/>

# 1-3-HamConnOrLaceable.gap

Read("LKH.gap");
Read("UndirectedGeneratingSets.gap");
Read("BipSetsOfCay.gap");
Read("IsHamiltonianCycle.gap");

IsHamPathToAllEndpoints := function(X, v, Endpoints)
	# returns true if X has a hamiltonian path from v 
	# to every vertex in Endpoints (other than v)
	local w, H;
	for w in Endpoints do
		if w <> v then
			# a hamiltonian path from v to w is the same as a hamiltonian cycle
			# that contains the edge [v,w] (in the graph that is obtained from X
			# by adding the edge vw).  As a signal to LKH that the edge vw is
			# undirected, we include its reverse in both the set of additional
			# edges and the set of required edges.
			H := LKH( X, [[v,w],[w,v]], [[v,w],[w,v]]);
			if 
				(H = fail)
					or 
				not(IsHamiltonianCycle( X, H, [[v,w],[w,v]], [[v,w],[w,v]] )) 
			then 
				return false;
			fi;
		fi;
	od;
	return true;
end;

# This function assumes that Cay(G;S) is bipartite, and that BipSetOfId is the bipartition
# set that contains the identity element. It returns a list of all the generating sets T 
# of G such that
# 	1. T can be obtained by adding a single nonidentity element g of G to S,
# 	2. Cay(G;T) is not bipartite (i.e., g is in BipSetOfId),
# 	3. g <= g^-1 (because we can replace an element by its inverse), 
# 	4. no proper subset of T gives a connected, non-bipartite Cayley graph of G.
# Actually, instead of returning all such generating sets, we choose one from each
# orbit of Aut(G).
NonBipGenSets := function (G, S, BipSetOfId)
	local NonBipElts, GenSets, g, T, Tless;
	NonBipElts := Filtered(BipSetOfId, g -> (g <> Identity(G)) and (g <= g^-1) );	
	GenSets := [];
	for g in NonBipElts do
		T := Concatenation(S, [g]);
		if
			ForAll(S, 
				function(s)
					# returns true if the subset Tless obtained by deleting s from T
					# does not give a connected, non-bipartite Cayley graph of G
					local Sless;
					Tless := Difference(T, [s]);
					return 
						Subgroup(G, Tless) <> G 
						or
						IsBipartite(CayleyGraph(G, Tless));
				end)
		then
			# add T to the list of generating sets, IF it is not already there
			if not(IsInListUpToUndirAut(T, GenSets)) then
				Add( GenSets, SortedList(T) );
			fi;
		fi;
	od;
	return GenSets;
end;
	
## This program checks whether all Cayley graphs of valency >= 3 are 
## hamiltonian connected or hamiltonian laceable (up to order 63).

# CayGS will be an abbreviation for CayleyGraph(G,S)
# initialize it so that GAP knows CayGS is a global variable
CayGS := 0; 

for k in [3..63] do
	Print(
		"\n\n----------------\n\n", 
		"k = ", k, ". ", 
		"There are ", NumberSmallGroups(k), " groups of this order.\n\n"
		);
	for GapId in [1..NumberSmallGroups(k)] do
		G := SmallGroup(k, GapId);
		Print("\n    G = SmallGroup(", k, ",", GapId, ") = ", StructureDescription(G), "\n");	
		IrredGenSets := IrredUndirGenSetsUpToAut(G);
		Print("        This group has ", Length(IrredGenSets), " irredundant generating sets.\n");
		
		# make a list of the irredundant generating sets whose Cayley graphs have valence 2
		ValenceTwoGenSets := Filtered(
			IrredGenSets, 
			S -> Length(AsSSortedList(Concatenation(S, List(S, s -> s^-1)))) = 2
			);
		Print("            ", Length(ValenceTwoGenSets), " of them have valence 2.\n");
		# remove these generating sets from consideration -- we know that their Cayley 
		# graphs are not hamiltonian connected/laceable
		GenSets := Difference(IrredGenSets, ValenceTwoGenSets);
		# for each generating set of valence 2, we construct the generating sets obtained 
		# by adding another element g of G (with g <= g^-1). We then add these generating 
		# sets to the list of those that need to be considered
		for S in ValenceTwoGenSets do
			for g in G do
				if (g <= g^-1) and not(g = Identity(G)) and not(g in S) and not(g^-1 in S) then
					Add(GenSets, Concatenation( S, [g] ) );
				fi;
			od;
		od;
		Print("        There are ", Length(GenSets), " generating sets to consider.\n");

		for S in GenSets do
			Print("          ", Position(GenSets, S), ". S = ", S, "\r");
			
			CayGS := CayleyGraph(G,S);
		
			# If Cay(G;S) is bipartite, let Endpoints be the bipartition set that does 
			# not contain Id.  Otherwise let Endpoints be all of G. We should be able to
			# find a hamiltonian path from Identity(G) to each nonidentity element of the
			# set Endpoints.
			BipSets := BipSetsOfCay(G,S);
			if BipSets = fail then # Cay(G;S) is not bipartite
				Endpoints := VertexNames(CayGS);
			else
				Endpoints := BipSets[2];
				BipSetOfId := BipSets[1]; # bipartition set that contains the identity element
			fi;				
			# We do not need to find paths to both x and the inverse of x
			Endpoints := Filtered(Endpoints, x -> (x <= x^-1));
		
			if IsHamPathToAllEndpoints(CayGS, Identity(G), Endpoints) then
				# Cay(G;S) is hamiltonian connected or hamiltonian laceable
				Print("          ", Position(GenSets, S), ". S = ", S, " is hamiltonian ");
				if BipSets = fail then
					Print("connected.\n");
				else
					Print("laceable\n");
				fi;
			else
				Print("\n");
				Error("\nUnable to confirm that this Cayley graph is hamiltonian connected/laceable\n");
			fi;
	
			# if Cay(G;S) is bipartite, it may become nonbipartite when another generator
			# is added to S. We check all possibilities T (up to automorphisms of G).
			if BipSets <> fail then # Cay(G;S) is bipartite
				TList := NonBipGenSets(G, S, BipSetOfId);
				Print("             with ", Length(TList), " nonbipartite extensions T:\n");
				for T in TList do;
					Print("               ", Position(TList, T), ". T = ", T, "\r");
					
					if IsHamPathToAllEndpoints( 
						CayleyGraph(G,T), 
						Identity(G), 
						Filtered(BipSetOfId, x -> x <= x^-1) 
						) 
					then
						Print("               ", Position(TList, T), ". T = ", T, " is hamiltonian connected.\n");
					else
						Print("\n");
						Error("Unable to confirm that this Cayley graph is hamiltonian connected\n");
					fi; 
				od;
			fi;
		od;
	od;
od;

Print("\n\nSuccess: all of the Cayley graphs are hamiltonian connected or hamiltonian laceable\n");

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

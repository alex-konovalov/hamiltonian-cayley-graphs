# SeveralHamCycs.gap

Read("HamiltonianCycles.gap");
Read("BipSetsOfCay.gap");

NumHamCycsToFind := function(GBar, SBar)
	# how many hamiltonian cycles should be returned by SeveralHamCycsInCay.
	# These numbers were determined empirically to have enough hamiltonian cycles 
	# for a nonzero voltage, without being excessively large

	local k, GapId;
	
	k := Order(GBar);
	if k = 40 then 
		GapId := IdSmallGroup(GBar)[2];
	fi;
	
	if k < 20 then
		return 20;
	elif k < 24 then
		return 30; # need 29 for D20
	elif k < 28 then
		return 50; # need 42 for C2 x C2 x S3
	elif k = 30 then
		return 130; # need 122 for D30
	elif (k = 40) and (GapId = 9) and (SBar = [GBar.1, GBar.3, GBar.1*GBar.2]) then
		return 4; # finding 4 ham cycles is fast, but the 5th takes a long time
	elif (k = 40) and (GapId = 9) then
		return 100;
	elif (k = 40) then
		return 25;
	elif (k = 42) then
		return 50; # need 46 for order 42
	elif (k = 44) then
		return 60; # need 53 for D44
	else
		return 40;
	fi;
end;	

Cay := function(G, S)
	# construct the (undirected) Cayley graph of the group G with the generating set S, but
	# multiplying on the RIGHT by elements of S. (The result is returned in GRAPE format.)
	local X, Spm;

	Spm := AsSSortedList( Concatenation( S, List(S, x -> x^(-1)) ));
    
    X := Graph( 
    	Group(()), 
    	[1..Order(G)], 
    	OnPoints, 
    	function(x,y) return (List(G)[x]^(-1) * List(G)[y] in Spm); end,
    	true
    	);    
    AssignVertexNames(X, List(G));   
    return X;   
end;

SeveralHamCycsInCay := function(GBar,SBar)
	# return a list of several hamiltonian cycles in CayleyGraph(GBar,SBar)
	# (The NumHamCycsToFind function prescribes the number of cycles.)

	local NumberOfHamCycsWanted, NumFound, HamCycs, NoRequiredEdges, 
	H, X, NeedToInvert, k, GapId;
	
	NumberOfHamCycsWanted := NumHamCycsToFind(GBar,SBar);
	NumFound := 0;
	HamCycs := [];
	NoRequiredEdges := [];
	k := Order(GBar);
	
	# there are a few Cayley graphs on which the HamiltonianCycles function is very slow 
	# but goes much faster if we use Cay(G,S) instead of CayleyGraph(G,S). But then the 
	# inversion operation needs to be applied, to turn right-multiplication into
	# left-multiplication, so we have a hamiltonian cycle in CayleyGraph(G,S) 
	if k in [40, 42, 44] then 
		GapId := IdSmallGroup(GBar)[2];
	fi;
	if 
		( (k = 40) and not(GapId in [1,2,3,12]) ) 
			or 
		( (k = 42) and (GapId in [3,4,5,6]) ) 
			or
		( (k = 44) and (GapId in [1]) and (SBar = [ GBar.1, GBar.3 ]) )
	then
		X := Cay(GBar, SBar);
		NeedToInvert := true;
	else
		X := CayleyGraph(GBar, SBar);
		NeedToInvert := false;
	fi;
	
	Print("           Found ", 0, " of ", NumberOfHamCycsWanted, " hamiltonian cycles\r");
	for H in HamiltonianCycles( X, NoRequiredEdges) do
		if NeedToInvert then
			H := List(H, x -> x^-1);
		fi;
		Add(HamCycs, H);
		NumFound := NumFound + 1;
		Print("           Found ", NumFound, " of ", NumberOfHamCycsWanted, " hamiltonian cycles\r");
		if NumFound >= NumberOfHamCycsWanted then
			break;
		fi;
	od;
	Print("\n");
		
	return HamCycs;
end;


SeveralHamCycsInRedundantCay := function(GBar, S0Bar, a)
	local s, a0, n, 
		AddEdges, g, VoltageOfH, HamCycs, H, i, Twists, GCDs, conductors,
		BadTwists, twist, numWanted, numDone, S;
	
	# For each abelian character, we ask LKH to find a hamiltonian cycle H that uses 
	# the edge [e,a], but no other edge whose endpoints are in the bipartition set of 
	# Cay(GBar,S0Bar) that contains the identity. Therefore, H has only one other a-edge. 
	# And we restrict this other a-edge to be one whose voltage will not cancel the voltage
	# of [e,a]. If LKH fails to find such a hamiltonian cycle, then we add this abelian
	# character to the list of bad twists that will have to be considered in the second
	# stage of the algorithm.
	HamCycs := [];
	BadTwists := [];
	for twist in AbelianCharacters(GBar) do
		AddEdges := [ [Identity(GBar), a] ];
		for g in BipSetsOfCay(GBar,S0Bar)[2] do
			if twist(g) <> -1 then
				Add(AddEdges, [g, a*g]);
			fi;
			if twist(g*a^(-1)) <> 1 then
				Add(AddEdges, [g, a^(-1)*g]);
			fi;
		od;
		H := LKH( CayleyGraph(GBar,S0Bar), AddEdges,[ AddEdges[1] ] );
		if (H <> fail) and IsHamiltonianCycle( CayleyGraph(GBar, Concatenation(S0Bar, [a])), H, [], []) then
			Add(HamCycs, H);
		else
			Add(BadTwists, StructuralCopy(twist));
		fi;
	od;
	
	# we are now in the second stage of the algorithm
	# We calculate the voltage for each hamiltonian cycle that we have found, and each 
	# bad twist. For each bad twist, we keep track of the GCD of the voltages for that
	# twist.
	numWanted := 100;
	numDone := 0;
	S := List(S0Bar, s -> [0,s]);
	Add(S, [1,a]);
	GCDs := List([1..Length(BadTwists)], x -> 0);
	conductors := List( BadTwists, twist -> Conductor( List( S0Bar, s -> twist(s) )));
	for H in HamCycs do
		for i in [1..Length(BadTwists)] do
			VoltageOfH := VoltageOfHamCycle(H,S,BadTwists[i]);
			if (VoltageOfH <> 0) then					
				GCDs[i] := Gcd(GCDs[i], Norm(VoltageOfH));
			fi;
		od;
	od;
	
	# if, for every bad twist, the GCD of the voltages has no prime factors bigger 
	# than 20 (that support the twist), then we have enough hamiltonian cycles
	if 
		ForAll([1..Length(BadTwists)], j -> 
			(GCDs[j] <> 0) 
				and
			(Maximum(Concatenation([1], Filtered(PrimeDivisors(GCDs[j]), p -> ( ((p-1) mod conductors[j]) = 0 )))) 
				<= 20
				)
			)
	then
		# done
	else
	
		if Length(S0Bar) = 1 and IsEvenInt(Order(GBar)) and (a in [S0Bar[1]^2, S0Bar[1]^-2]) then
			# we constructed a hamiltonian cycle by hand for this special case
			s := S0Bar[1];
			a0 := s^2;
			n := Order(GBar);
			H := Concatenation( List([0..n/2-1], i -> a0^i), [s^-1], List([1..n/2-1], i -> s^-1*a0^-i), [Identity(GBar)]);
			Add(HamCycs, H);
			for i in [1..Length(BadTwists)] do
				VoltageOfH := VoltageOfHamCycle(H,S,BadTwists[i]);
				if (VoltageOfH <> 0) then					
					GCDs[i] := Gcd(GCDs[i], Norm(VoltageOfH));
				fi;
			od;
		fi;
		
		# try again
		GCDs := List([1..Length(BadTwists)], x -> 0);
		conductors := List( BadTwists, twist -> Conductor( List( S0Bar, s -> twist(s) )));
		for H in HamCycs do
			for i in [1..Length(BadTwists)] do
				VoltageOfH := VoltageOfHamCycle(H,S,BadTwists[i]);
				if (VoltageOfH <> 0) then					
					GCDs[i] := Gcd(GCDs[i], Norm(VoltageOfH));
				fi;
			od;
		od;
		if 
			ForAll([1..Length(BadTwists)], j -> 
				(GCDs[j] <> 0) 
					and
				(Maximum(Concatenation([1], Filtered(PrimeDivisors(GCDs[j]), p -> ( ((p-1) mod conductors[j]) = 0 )))) 
					<= 20
					)
				)
		then
			# done
		else
			# look for additional hamiltonian cycles until all prime divisors are small
			for H in HamiltonianCycles( CayleyGraph(GBar,Concatenation(S0Bar, [a])), [[Identity(GBar),a]]) do
				Add(HamCycs, H);
				for i in [1..Length(BadTwists)] do
					VoltageOfH := VoltageOfHamCycle(H,S,BadTwists[i]);
					if (VoltageOfH <> 0) then					
						GCDs[i] := Gcd(GCDs[i], Norm(VoltageOfH));
					fi;
				od;
				if ForAll([1..Length(BadTwists)], j -> (GCDs[j] <> 0) and 
					(Maximum(Concatenation([1], Filtered(PrimeDivisors(GCDs[j]), p -> ( ((p-1) mod conductors[j]) = 0 )))) <= 5))
				then
					break;
				fi;
				numDone := numDone + 1;
				if numDone >= numWanted then
					break;
				fi;
			od;
		fi;
	fi;
	return HamCycs;
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

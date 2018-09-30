# CallLKHOnLiftsOfSBar.gap

# if we are unable to find a hamiltonian cycle in CayleyGraph(GBar,SBar) that has 
# nonzero voltage modulo p, then we call LKH on all possible Cayley graphs that 
# have CayleyGraph(GBar,SBar) as their quotient. That is, we construct all possible
# lifts S of SBar to G (up to conjugacy), and call LKH on each Cayley graph.

Read("LKH.gap");

MakeZpxGBar := function(GBar, twist, p)
	# construct the semidirect product of GBar with the cyclic group of order p
	# (with respect to the twist function obtained by converting roots of unity
	# to the corresponding automorphisms of the cyclic group). The result is
	# returned as a GAP group.
	# The function returns a pair consisting of:
	#	* the semidirect product G
	#	* a function that translates ordered pairs [z,gbar] to elements of G
	local Zp, AutZp, gens, ConductorOfTwist, base, twistp, G, PairToG;
	Zp := CyclicGroup(p);
	AutZp := AutomorphismGroup(Zp);
	gens := GeneratorsOfGroup(GBar);
	
	ConductorOfTwist := Conductor( List( GBar, x -> twist(x) ) );
	if not( ((p-1) mod ConductorOfTwist) = 0 ) then
		# the prime p does not support this twist
		return fail;
	elif IsOddInt(ConductorOfTwist) then
		# need to allow E(2) = -1 even if the conductor is odd
		ConductorOfTwist := 2 * ConductorOfTwist;
	fi;
	
	# translate twist (whose values are roots of unity) to a function
	# with values in Aut(Zp), so we can form the semidirect product
	base := MinimalGeneratingSet(AutZp)[1]^((p-1)/ConductorOfTwist);
	twistp := function(g)
		local descrip;
		descrip := DescriptionOfRootOfUnity( twist(g) );
		return (base^(ConductorOfTwist/descrip[1]))^(descrip[2]);
	end;
	 
	G := SemidirectProduct(
		GBar, 
		GroupHomomorphismByImages(GBar, AutZp, gens, List(gens, x -> twistp(x))), 
		Zp);
	
	PairToG := function(x)
			local iG, iZ;
			iG := Embedding(G, 1); # iG is the embedding of GBar into G
			iZ := Embedding(G, 2); # iZ is the embedding of Zp into G
			return ( Image(iZ, (Zp.1)^x[1]) * Image(iG,x[2]) );
		end;
	return [G, PairToG];
end;	


CallLKHOnLiftsOfSBar := function( GBar, SBar, twist, p )
	# ask LKH to find hamiltonian cycles in all possible lifts of Cay(GBar,SBar)
	local Zp, AutZp, j, base, twistp, gens, G, iG, iZ,
		toG, z, h, H, SList, S, ZpxGBar, HowMuchOfZp, CayGS;
	
	if IsEvenInt(p) or not(IsPrime(p)) then
		Error("p must be an odd prime, not ", p, "\n");
	fi;
	
	# construct the semidirect product G that covers GBar
	ZpxGBar := MakeZpxGBar(GBar, twist, p); 
	if ZpxGBar = fail then
		Print("The prime ", p, " does not support this twist.\n");
		return true;
	fi;
	G := ZpxGBar[1]; # G is now the semidirect product of CyclicGroup(p) with GBar
	toG := ZpxGBar[2];
	
	HowMuchOfZp := function(SBar, twist, i)
	# To lift an element of GBar to an element of G, we multiply it by an element of Zp.
	# However, elements of SBar usually do not need to be multiplied by every element of Zp.
	# For each i, this function returns that subset of Zp that SBar[i] needs to be 
	# multiplied by.
	
		if ( (Order(SBar[i]) = 2) and (twist(SBar[i]) = 1) ) 
		then # Factor Group Lemma applies if SBar[i]  projects nontrivially to Zp
			return List( [ [0, SBar[i] ] ], toG);
		
		elif ( (twist(SBar[i]) <> 1)  and ForAll([1..i-1], j -> (twist(SBar[i]) = 1)))
		then # conjugate by an element of Zp to kill the first element with a nontrivial twist
			return List( [ [0, SBar[i] ] ], toG);

		elif
			ForAll([1..i-1], j ->
				( ( (Order(SBar[j]) = 2) and (twist(SBar[j]) = 1) )
			or
			( (twist(SBar[j]) <> 1)  and ForAll([1..j-1], k -> (twist(SBar[k]) = 1) ))))
		then # we can normalize the first nonzero projection to Zp to be 1
			return List( Cartesian([ [0, 1], [SBar[i]] ]), toG );

		else # we need to allow all elements of Zp
			return List( Cartesian([ [0..p-1], [SBar[i]]]), toG );
			
		fi;
	end;
	
	# Each lift S is constructed by choosing an element of HowMuchOfZp(i) for each i.
	# Thus, the set of all possible S is the Cartesian product
	SList := Cartesian( List( [1..Length(SBar)], i -> HowMuchOfZp(SBar, twist, i)) );
	# However, we are only interested in sets that generate G:
	SList := Filtered( SList, S -> (Subgroup(G,S) = G) );
	Print(
		"                        There are ", Length(SList), 
		" lifts of SBar to G = ", StructureDescription(G), "\n"
		);

	for S in SList do
		Print("                            ", Position(SList, S), ". S = ", S, "\r");
		CayGS := CayleyGraph(G, S);
		H := LKH(CayGS, [], []);
		if H <> fail and IsHamiltonianCycle(CayGS, H, [], []) then
			Print(
				"                            ", Position(SList, S), ". S = ", S, 
				": LKH found a hamiltonian cycle.\n"
			);
		else
			Error("\nLKH did not find a hamiltonian cycle: ", S, "\n");
		fi;
	od;
	return true;
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

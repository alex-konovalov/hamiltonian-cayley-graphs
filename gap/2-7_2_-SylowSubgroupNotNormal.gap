# 2-7(2)-SylowSubgroupNotNormal.gap

Read("LKH.gap");
Read("UndirectedGeneratingSets.gap");

for k in [1..47] do
	LargestPrimeDivisor := Maximum(Concatenation(PrimeDivisors(k), [1]));
	Print("\nk = ",k, ": largest prime divisor is ", LargestPrimeDivisor, ".\n");
	for p in [LargestPrimeDivisor+1..k] do
		if (IsPrime(p)) and ForAny(DivisorsInt(k), d -> (d > 1) and ((d mod p) = 1) ) then
			Print("\n    p = ", p, ": Checking the ", NumberSmallGroups(k*p)," groups of order ", k, " x ", p, " = ", k*p, ".\n");
			for GapId in [1..NumberSmallGroups(k*p)] do
				G := SmallGroup(k*p,GapId);
				P := SylowSubgroup(G,p);
				if not(IsNormal(G,P)) then
					Print("        G = SmallGroup(", k*p, ",", GapId, ") = ", StructureDescription(G), "\r"); 
					GenSets := IrredUndirGenSetsUpToAut(G);
					Print(
						"        G = SmallGroup(", k*p, ",", GapId, ") = ", StructureDescription(G), 
						" has ", Length(GenSets), " generating sets.\n"
						);
					for S in GenSets do
						Print("                ", Position(GenSets, S), ". S = ", S, ": \r");
						H := LKH( CayleyGraph(G,S), [], []);
						if 
							H <> fail 
								and
							IsHamiltonianCycle( CayleyGraph(G,S), H, [], [] )
						then
							Print(
								"                ", Position(GenSets, S), ". S = ", S, ": ",
								" LKH found a hamiltonian cycle\n"
								);
						else
							Print("\n");
							Error(" LKH failed to find a hamiltonian cycle!\n");
						fi;
					od;
				fi;
			od;
		fi;
	od;
od;
Print("\nSuccess: LKH found a hamiltonian cycle whenever Zp is not normal in G.\n");

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

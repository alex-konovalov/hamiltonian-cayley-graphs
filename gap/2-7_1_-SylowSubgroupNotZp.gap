# 2-7(1)-SylowSubgroupNotZp.gap

Read("LKH.gap");
Read("UndirectedGeneratingSets.gap");

# keep a list of the orders of groups that arise
# (for printing in the summary at the end)
GroupOrders := [];

for k in [1..47] do
	p := Maximum(Concatenation(PrimeDivisors(k), [1]));
	Print("\nk = ",k,": The largest prime divisor of k is ", p,".\n");
	if k*p = 64 then
		Print("        Connected Cayley graphs of order 64 are hamiltonian.\n");
	else
		NumberOfGroups := NumberSmallGroups(k*p);
		Print(
			"    There are ", NumberOfGroups," groups of order ",
			 k, " x ", p," = ", k*p, ".\n"
			);
		for GapId in [1..NumberOfGroups] do
			G := SmallGroup(k*p,GapId);
			if IsAbelian(G) then
				Print(
					"\n        ", GapId, ". SmallGroup(", k*p, ",", GapId, ") = ",
					StructureDescription(G), " is abelian\n"
					);
			elif (k*p = 1058) and (GapId = 4) then
				# special case where IrredUndirGenSetsUpToAut may run out of memory
				Print(
					"\n        ", GapId, ". SmallGroup(", k*p, ",", GapId, ") = ",
					StructureDescription(G), " is dihedral type of order 2p^2\n"
					);
			else
				Print(
					"\n        ", GapId, ". SmallGroup(", k*p, ",", GapId, ") = ",
					StructureDescription(G), "\r"
					);
				GenSets := IrredUndirGenSetsUpToAut(G);
				Print(
					"        ", GapId, ". SmallGroup(", k*p, ",", GapId, ") = ",
					StructureDescription(G), 
					" has ", Length(GenSets), " generating sets.\n"
					);
				for S in GenSets do
					Print("            ", Position(GenSets,S), ". S = ", S, ": \r");
					CayGS := CayleyGraph(G,S);
					H := LKH(CayGS, [],[]);
					if 
						(H <> false) 
							and 
						IsHamiltonianCycle(CayGS, H, [], []) 
					then
						Print(
							"            ", Position(GenSets,S), ". S = ", S,
							": LKH found a hamiltonian cycle\n"
							);
					else
						Print("\n");
						Error("LKH failed to find a hamiltonian cycle!\n");
					fi;
				od;
			fi;
		od;
	fi;
	Add(GroupOrders,[k, k*p]);
od;
Print("\nSuccess: All connected Cayley graphs of the following orders are hamiltonian:\n");
for kPair in GroupOrders do
	Print("    ", kPair[2], " (k = ", kPair[1], ")\n");
od;

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

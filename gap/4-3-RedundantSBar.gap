# 4-3-RedundantSBar.gap

Read("UndirectedGeneratingSets.gap");
Read("IsHamiltonianCycle.gap");
Read("BipSetsOfCay.gap");
Read("AbelianCharacters.gap");
Read("VoltageOfHamCycle.gap");
Read("CallLKHOnLiftsOfSBar.gap");

# do not need to verify:
Read("SeveralHamCycs.gap");


# initialize these variables so GAP knows that they are global
GBar := [];
ConductorOfTwist := 0;

for k in [1..47] do
	if IsOddInt(k) then
		Print("\n\nk = ", k, ". No nontrivial bipartite Cayley graphs have odd order.\n");
	else
		Print("\n\nk = ", k, ". There are ", NumberSmallGroups(k)," groups of order ", k, ".\n");
		for GapId in [1..NumberSmallGroups(k)] do
			GBar := SmallGroup(k,GapId);
			Print("\n    GBar = SmallGroup(", k, ", ", GapId, ") = ", StructureDescription(GBar), "\r");

			# loop through irredundant generating sets of GBar
			GenSets := IrredUndirGenSetsUpToAut(GBar);
			Print(
				"    GBar = SmallGroup(", k, ", ", GapId, ") = ", StructureDescription(GBar), 
				" has ", Length(GenSets), " irredundant generating sets\n"
				);
			for S0Bar in GenSets do 
				Print("        ", Position(GenSets, S0Bar), ". S0Bar = ", S0Bar, ":\r");
			
				BipSets := BipSetsOfCay(GBar,S0Bar);
				if BipSets = fail then
					Print("        ", Position(GenSets, S0Bar), ". S0Bar = ", S0Bar, ": not bipartite\n");
				else

					# Cay(GBar, S0Bar) is bipartite, so we need to look at the nonbipartite extensions.
					# The possibilities for a are the nonidentity elements of the bipartition set that
					# contains the identity, but we do not need to consider both a and a^-1.
					SetOfa := Filtered(Difference( BipSets[1], [Identity(GBar)] ), x -> (x <= x^-1) );
				
					for a in SetOfa do

						Print(
							"        ", Position(GenSets, S0Bar), ". S0Bar = ", S0Bar, ": trying ", 
							Position(SetOfa, a), " of ", Length(SetOfa), " (a = ", a, ")\r"
							);
					
						# construct the corresponding generating set S of G			
						SBar := Concatenation( S0Bar, [a]);
						S := Concatenation( List(S0Bar, s -> [0,s]), [ [1,a] ]);
															
						# make a list of a few hamiltonian cycles in Cay(GBar, SBar)
						HamCycs := SeveralHamCycsInRedundantCay(GBar, S0Bar, a);
					
						# loop through all abelian characters of GBar
						for twist in AbelianCharacters(GBar) do
							ConductorOfTwist := Conductor( List( GBar, x -> twist(x) ) );
						
							if (Order(a) = 2) and (twist(a) = 1) then
								# Lemma 2.13(9) rules out this case, so we ignore it
							else
								# Verify that none of the elements of S are redundant.
								# (Note that if a is the unique element of SBar \ {t} with a 
								# nontrivial twist, then S \ {t} does not generate G, even if 
								# it generates GBar.)
								if 
									ForAll(S0Bar, t -> 
										not(Subgroup(GBar, Difference(S0Bar, [t]) ) = GBar)
											or 
										( 
											(twist(a) <> 1) 
												and 
											ForAll( Difference(S0Bar, [t]), s -> (twist(s) = 1) )
										)
									)
								then
									# we calculate the GCD of (the norms of) the voltages in our 
									# list of hamiltonian cycles (after verifying that they are 
									# indeed hamiltonian cycles)
									GCD := 0;
									for H in HamCycs do
										if IsHamiltonianCycle( CayleyGraph(GBar,SBar), H, [], []) then
											VoltageOfH := VoltageOfHamCycle(H, S, twist);
											if VoltageOfH <> 0 then
												GCD := Gcd(GCD, Norm(VoltageOfH));
											fi;
										fi;
									od;
								
									# if p is not a divisor of the GCD, then there is some
									# hamiltonian cycle whose voltage is nonzero mod p. So
									# Cay(G,S) has a hamiltonian cycle (cf. Lemma 3.3 with 
									# n = 1). Thus, we now only need to consider the primes
									# that divide the GCD. Furthermore, we only need to 
									# consider primes that support the twist that we are 
									# looking at now (and are larger than the largest prime 
									# factor of k).
									
									if GCD = 0 then
										Error("GCD of voltages is 0");
										
									elif Maximum(Concatenation([1], PrimeDivisors(GCD))) > 10^8 then
										Error(
											"The prime p = ", Maximum(PrimeDivisors(GCD)), 
											"is far too large for LKH to handle (and GAP",
											" may not guarantee that it is prime)."
											);
											
									elif
										ForAll(PrimeDivisors(GCD), p-> 
											(
												not( ((p-1) mod ConductorOfTwist) = 0 )
													or 
												p <= Maximum(PrimeDivisors(Order(GBar))) 
											) 
										) 
									then
									
										# GCD is not divisible by any primes that support this twist
										# so we can move on to the next twist
									else	
										# some prime divisors of the GCD support this twist 
										Print(
											"\n            GCD of voltages is ", GCD, 
											" for twist(S0Bar) = ", List(S0Bar, twist),
											" and twist(a) = ", twist(a), ".\n",
											"            Prime divisors are ", PrimeDivisors(GCD), "\n"
											);
										for p in PrimeDivisors(GCD) do
									
											if p <= Maximum(PrimeDivisors(k)) then
												Print("                p = ", p, " is not larger than",
													" the largest prime divisor ", Maximum(PrimeDivisors(k)), 
													" of k = ", k, "\n"
													);
												
											elif not( ((p-1) mod ConductorOfTwist) = 0 ) then
												Print(
													"                p = ", p, " does not support",
													" this twist of conductor ", ConductorOfTwist, "\n"
													);
												
											else
												Print("                p = ", p, ": we call LKH\n");
												CallLKHOnLiftsOfSBar( GBar, SBar, twist, p );
											fi;
										od;	
									fi;
								fi;
							fi;
						od;
					od;
					Print("        ", Position(GenSets, S0Bar), ". S0Bar = ", S0Bar, ": all ", Length(SetOfa), " choices of a are ok", "          \n");
				fi;
			od;
		od;
	fi;
	Print("\n--------------------------------------------------------------------\n");
	Print("--------------------------------------------------------------------\n");
od;

Print(
	"Success: found hamiltonian cycles in Cay(G,S) whenever SBar is redundant",
	" and not bipartite, but S0Bar is bipartite"
	);

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

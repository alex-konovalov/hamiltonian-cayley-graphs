# 3-4-IrredundantSBar.gap

Read("UndirectedGeneratingSets.gap");
Read("AbelianCharacters.gap");
Read("BasisGeneratingSets.gap");
Read("CheckForIrredundantSpecialCase.gap");
Read("VoltageOfHamCycle.gap");
Read("CallLKHOnLiftsOfSBar.gap");
Read("CyclotomicToFp.gap");
Read("Det.gap");

# do not need to verify:
Read("SeveralHamCycs.gap");
Read("FindNonzeroDet.gap");


VoltageRows := []; # initialize so GAP knows it is a global variable

for k in [1..47] do
	Print(
		"\nk = ", k, ": Largest prime divisor is ", Maximum(Concatenation([1], PrimeDivisors(k))), 
		". There are ", NumberSmallGroups(k), " groups of this order.\n"
		);

	for GapId in [1..NumberSmallGroups(k)] do 
		GBar := SmallGroup(k,GapId);
		Print("\n    GBar = SmallGroup(", k, ", ", GapId, ") = ", StructureDescription(GBar), "\r");

		# loop through the irredundant generating sets of GBar
		GenSets := IrredUndirGenSetsUpToAut(GBar);
		Print(
			"    GBar = SmallGroup(", k, ", ", GapId, ") = ", StructureDescription(GBar), 
			" has ", Length(GenSets), " irredundant generating sets\n"
			);
		for SBar in GenSets do 
			Print(
				"\n        ", Position(GenSets, SBar), ". SBar = ", SBar, 
				"   in GBar = SmallGroup(", k, ", ", GapId, ") = ", StructureDescription(GBar), "\n"
				);
			
			# make a list of a few hamiltonian cycles in CayleyGraph(GBar,SBar)
			SeveralHamCycsInCG := SeveralHamCycsInCay(GBar,SBar);
			
			# loop through all abelian characters of GBar
			Twists := AbelianCharacters(GBar);
			Print("           There are ", Length(Twists), " twists to form semidirect products.\n");
			for twist in Twists do
				Print("\n            ", Position(Twists, twist), ". Twist = ", List(SBar, s -> twist(s)), "\r");
				ConductorOfTwist := Conductor( List( GBar, x -> twist(x) ) );
										
				# check whether Lemma 3.1 (or other theoretical result) applies
				IrredundantSpecialCase := CheckForIrredundantSpecialCase(GBar, SBar, twist, k, GapId);
				if IrredundantSpecialCase <> fail then
					# some lemma provides a hamiltonian cycle, so do not need to do anything else
					Print(
						"            ", Position(Twists, twist), ". Twist = ", List(SBar, s -> twist(s)),
						": ", IrredundantSpecialCase, "\n"
						);

				else 
					# no special case applies, so we try to find hamiltonian cycles 
					# whose voltages make a nonzero determinant that will let us 
					# apply Lemma 3.3
				
					Print("\n");
					
					# let SList be the list of generating sets S_1,S_2,... specified in Lemma 3.3
					SList := BasisGeneratingSets(SBar, twist);

					# VoltageRows will be a list of the rows that are candidates to be put into 
					# the matrix whose determinant appears in the statement of Lemma 3.3
					VoltageRows := []; 
					
					for H in SeveralHamCycsInCG do
						# confirm that H really is a hamiltonian cycle
						if not(IsHamiltonianCycle( CayleyGraph(GBar, SBar), H, [], [] ) ) then
							Error("\nNot a hamiltonian cycle: ", H, "\n\n");
						fi;
						
						# the hamiltonian cycle H has several associated voltages (corresponding
						# to the elements of SList). These voltages can be used as one row of 
						# the matrix whose determinant must be nonzero to apply Lemma 3.3.
						OneRowOfVoltages := List(
							SList,
							S -> VoltageOfHamCycle(H,S,twist) 
							);
						# add this to the list of these possible rows for the determinant
						Add(VoltageRows,OneRowOfVoltages);
					od;
			
					HamCycsWithNonzeroDet := FindNonzeroDet(VoltageRows);
					if HamCycsWithNonzeroDet = fail then
						Error("unable to find a nonzero determinant\n");
					else
						
						NonzeroNormDet := Norm(Det( List(HamCycsWithNonzeroDet, i -> VoltageRows[i]) ));
						# verify that we do have a nonzero determinant
						if NonzeroNormDet = 0 then
							Error("Norm of determinant is unexpectedly zero\n");
						fi;
						Print(
							"                For hamiltonian cycles ", HamCycsWithNonzeroDet, 
							", norm(determinant) = ", NonzeroNormDet, "\n                ",
							"Prime divisors are ", PrimeDivisors(NonzeroNormDet), "\n"
							);
						
						# Lemma 3.3 provides a hamiltonian cycle for every prime p that does
						# NOT divide NonzeroNormDet. We now consider the remaining primes
						for p in PrimeDivisors(NonzeroNormDet) do
						
							if p > 10^8 then
								Error(
									"The prime p = ", p, "is far too large for LKH to handle",
									" (and GAP may not guarantee that it is prime)."
									);
					
							elif p <= Maximum( Concatenation( PrimeDivisors(k), [1] ) ) then
								Print("                    p = ", p, " is not larger than",
									" the largest prime divisor ", Maximum( Concatenation( PrimeDivisors(k), [1] )), 
									" of k = ", k, "\n"
									);
							
							elif ((p-1) mod ConductorOfTwist) <> 0 then
								Print(
									"                    p = ", p, " does not support",
									" this twist of conductor ", ConductorOfTwist, "\n"
									);
							
							else
								HamCycsWithNonzeroDetModp := FindNonzeroDetModp(VoltageRows, p);
								if HamCycsWithNonzeroDetModp <> fail then 
									DetModp := Det( 
										List(HamCycsWithNonzeroDetModp, 
											i -> List(VoltageRows[i], 
												voltage -> CyclotomicToFp(voltage, p))
											)
										);
									if DetModp = 0 then
										Error("Determinant is unexpectedly divisible by p\n");
									else
										Print(
											"                    p = ", p, ": for ",
											"hamiltonian cycles ", HamCycsWithNonzeroDetModp, 
											", determinant = ", DetModp, " is nonzero in Fp\n"
											);
									fi;
								else
									Print(
										"                    p = ", p, ": determinants ",
										"are 0 modulo p, so we call LKH\n"
										);
									CallLKHOnLiftsOfSBar( GBar, SBar, twist, p );
								fi;													
							fi;
						od;
					fi;				
				fi;
			od;
		od;
	od;
	Print("\n--------------------------------------------------------------------\n");
	Print("--------------------------------------------------------------------\n");
od;

Print("Success: found hamiltonian cycles in Cay(G,S) whenever SBar is irredundant");

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

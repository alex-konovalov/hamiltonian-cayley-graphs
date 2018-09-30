# 4-5-Valence2.gap

Read("UndirectedGeneratingSets.gap");
Read("AbelianCharacters.gap");
Read("VoltageOfHamCycle.gap");
Read("CallLKHOnLiftsOfSBar.gap");
Read("IsHamiltonianCycle.gap");

# do not need to verify:
Read("HamiltonianCycles.gap");

# initialize these variables so GAP knows that they are global
GBar := [];
ConductorOfTwist := 0;
S0Bar := 0;
twist := [];

for k in [4..47] do
	Print("\nk = ", k, ".\n");
	
	for CyclicOrDihedral in ["cyclic", "dihedral"] do
		if CyclicOrDihedral = "cyclic" then
			GBar := CyclicGroup(k);
			S0Bar := [GBar.1]; # the 1-element generating set of GBar
		elif IsEvenInt(k) then
			GBar :=  DihedralGroup(k);
			S0Bar := [GBar.1,GBar.1*GBar.2]; # two elements of order 2 that generate GBar
		else
			break; # there are no dihedral groups of odd order
		fi;

		Print("\n    GBar = ", StructureDescription(GBar)," = SmallGroup(",k,", ",IdSmallGroup(GBar)[2],"), S0Bar = ", S0Bar, "\n");

		if
			(Length(AsSSortedList( Concatenation(S0Bar, List(S0Bar, s -> s^-1)))) <> 2) # valence is not 2
				or
			(Subgroup(GBar, S0Bar) <> GBar) # S0Bar does not generate GBar
		then
			Error("This is the wrong generating set!\n"); 
		fi;
		
		SetOfa := Filtered(GBar,
			a -> 
				not(a = Identity(GBar)) 
					and 
				(not(a in S0Bar)) 
					and 
				(not(a^-1 in S0Bar)) 
					and 
				(a <= a^-1) 
			);
		for a in SetOfa do
			Print(
				"        trying ", Position(SetOfa, a), " of ", Length(SetOfa), 
				" (a = ", a, ")     \r"
				);
			SBar := Concatenation(S0Bar,[a]);

			# a is the only element of S with nonzero projection to Z
			S := List(SBar,sbar->[0,sbar]);
			Add(S,[1,a]);

			# make a list of a few hamiltonian cycles in Cay(GBar,SBar)
			HamCycs := [];
			NumberOfHamCycsWanted := 20;
#			Print("        found ", Length(SomeHamCycs), " of ", NumberOfHamCycsWanted, "hamiltonian cycles\r");
			for H in HamiltonianCycles( CayleyGraph(GBar,SBar), [[Identity(GBar),a]]) do
				
				Add(HamCycs, H);
#				Print("        found ", Length(HamCycs), " of ", NumberOfHamCycsWanted, "hamiltonian cycles\r");
				
				if Length(HamCycs) >= NumberOfHamCycsWanted then
					break;
				fi;
			od;
			
			for twist in AbelianCharacters(GBar) do
				if ( (Order(a) = 2) and (twist(a) = 1) ) then
					# this contradicts Lemma 2.13(9), so we do not have to do anything
				else
					ConductorOfTwist := Conductor( List( GBar, x -> twist(x) ) );
				
					# this is the same as in Proposition 4.3
					
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
			od;
					
		od;
		Print(
			"        trying ", Position(SetOfa, a), " of ", Length(SetOfa), 
			": all extensions of S0Bar are ok   \n"
			);

	od;
	Print("\n-------------------------------------------------------------------- \n");
	Print("-------------------------------------------------------------------- \n");
od;

Print(
	"Success: found hamiltonian cycles in Cay(G,S) whenever SBar is a redundant",
	" extension of generating set of valence 2"
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

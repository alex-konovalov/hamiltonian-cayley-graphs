# VoltageOfHamCycle.gap

# The function VoltageOfHamCycle computes the voltage in Gtilde of a hamiltonian cycle
# in a Cayley graph CayleyGraph(GBar,SBar). It is assumed that Gtilde is the semidirect
# product Z x GBar of GBar with a subring of the complex numbers, with respect to the 
# given twist function.  Elements of Gtilde are represented as ordered pairs [z, gbar].
#	S = a generating set of Gtilde
#  	H = a hamilitonian cycle in Cay(GBar,SBar), written as a list of elements of GBar

# computes the inverse of a given ordered pair x in Gtilde
InverseSemiDirect := function(x,twist) 
	return [ (-1) * twist(Inverse(x[2])) * x[1], Inverse(x[2]) ];
end;

# computes the product of two ordered pairs in Gtilde
ProductSemiDirect := function(x,y,twist)
	return [ x[1] + twist(x[2]) * y[1], x[2] * y[2] ];
end;

VoltageOfHamCycle := function(HBar, S, twist)
	# compute the voltage of the hamiltonian cycle HBar in CayleyGraph(GBar,SBar)
	local HBarEdges, HEdges, WhereWeAre, j, i, IdentityOfGBar;
	
	HBarEdges := [];
	HEdges := [];

	# Translate HBar from a list of elements of GBar to a list of elements of SBar
	for j in [1..Length(HBar)-1] do
		HBarEdges[j] := HBar[j+1] * HBar[j]^-1;
	od;

	# Translate HBarEdges from a list of elements of SBar to a list of elements of S 
	for j in [1..Length(HBarEdges)] do
		for i in [1..Length(S)] do
			if HBarEdges[j] = S[i][2] then
				HEdges[j] := S[i];
		    elif HBarEdges[j] = Inverse(S[i][2]) then	
				HEdges[j] := InverseSemiDirect(S[i],twist);
		    fi;
		od;
	od;
	
	# we will need the identity element of GBar
	# (empty list is not a ham cyc, so HBar[1] should always exist)
	IdentityOfGBar := HBar[1]*HBar[1]^-1; 
	
	# walk starts at the identity element of Gtilde
	WhereWeAre:= [ 0, IdentityOfGBar ];
	# walk through Gtilde, by using the successive edges in HEdges
	for i in [1..Length(HEdges)] do
		WhereWeAre:= ProductSemiDirect( HEdges[i], WhereWeAre, twist ); 
	od;
	# voltage is the Z component of the endpoint of the walk
	return WhereWeAre[1];
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


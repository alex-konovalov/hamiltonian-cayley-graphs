# BasisGeneratingSets.gap

BasisGeneratingSets := function(SBar, twist)
	# construct the list of generating sets S_1,S_2,... specified in Lemma 3.3
	local SList, B0IsNonempty, i, S;
	
	SList := [];
	B0IsNonempty := false; # we have not yet put anything into the set B_0
	
	# look at each element s of SBar. If s is one of the b_i's, then construct the
	# corresponding generating set containing [1,s], and add it to SList.
	# Otherwise, just ignore it, because either s is an a_i or s is in B_0
	for i in [1..Length(SBar)] do  
		# construct the generating set S in which only the ith entry has nonzero Zp-coordinate
		S := List(SBar, s->[0,s]);
		S[i] := [1, SBar[i]];
		
		if Order(SBar[i])=2 and twist(SBar[i]) = 1 then
			# SBar[i] is an a_i, so do not add S to the list
		elif twist(SBar[i]) <> 1 and B0IsNonempty = false then
			# this is the first element with a nontrivial twist, so we can put it in B_0
			# (and record the fact that we have already put an element in B_0)
			# and we do not add S to the list
			B0IsNonempty := true;
		else
			# s is one of the b_i's, so add S to the list
			Add(SList, S) ;
		fi;
	od;
	
	return SList;
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

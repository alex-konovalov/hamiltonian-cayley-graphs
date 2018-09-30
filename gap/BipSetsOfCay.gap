# BipSetsOfCay.gap

LoadPackage("grape");

BipSetsOfCay := function(G, S)
	# If Cay(G,S) is bipartite, then this function returns a list of 2 elements. The first
	# element of the list is the bipartition set that contains the identity element of G, 
	# the other element of the list is the other bipartition set. Elements of the two sets
	# are specified by their vertex names. If Cay(G,S) is not bipartite, then the function 
	# returns fail.
	local CayGS, Bicomps, BipSets;
	
	CayGS := CayleyGraph(G,S);
	Bicomps := Bicomponents( CayGS ); 

	if IsEmpty(Bicomps) then
		return fail;
	else # Cay(G,S) is bipartite
	
		# let BipSets[1] be the bipartition set that contains the identity element, and
		# let BipSets[2] be the other one.  Also convert each of them to a list of 
		# elements of GBar
		if Identity(G) in List(Bicomps[1], i -> VertexName(CayGS,i)) then
			BipSets := [
				List(Bicomps[1], i -> VertexName(CayGS,i)),
				List(Bicomps[2], i -> VertexName(CayGS,i))
				];
		else 
			BipSets := [
				List(Bicomps[2], i -> VertexName(CayGS,i)),
				List(Bicomps[1], i -> VertexName(CayGS,i))
				];
		fi;
		return BipSets;
	fi;
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

# CyclotomicToFp.gap
 
CyclotomicToFp := function(x, p)
	# this is the unique homomorphism from the ring of cyclotomics of conductor dividing p-1
	# to the finite field Fp, that takes E(p-1) to Z(p,1).

	# Letting N = p-1 (and cyc = x), the GAP manual entry 18.1.10 for CoeffsCyc tells us that
	#	x = CoeffsCyc(x, p-1) * List( [1..p-1], j -> E(p-1)^(j-1) );
	# We simply replace E(p-1) with Z(p,1) in this formula.
	
	return CoeffsCyc(x, p-1) * List( [1..p-1], j -> Z(p,1)^(j-1) );
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

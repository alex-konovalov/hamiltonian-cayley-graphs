# IsHamiltonianCycle.gap

LoadPackage("grape");

IsHamiltonianCycle := function(X, H, AdditionalEdges, RequiredEdges)
	# returns true if H is a hamiltonian cycle in the graph obtained from X
        # by adding the directed edges in the list AdditionalEdges, and H uses every
        # directed edge in the list RequiredEdges
        # (an undirected edge is specified by listing both it and its inverse)
    # H, elements of AdditionalEdges, and elements of RequiredEdges are given as lists
    # of names of vertices of X
    local EdgesOfX;
    
    # convert vertex numbers to vertex names in a list of the (directed) edges of X
    EdgesOfX := List(DirectedEdges(X), e -> [VertexName(X, e[1]), VertexName(X, e[2])]);
    
    if Length(AsSSortedList(VertexNames(X))) <> Length(VertexNames(X)) then
    	Error("IsHamiltonianCycle requires all vertex labels to be distinct.\n");
    	return false;
    fi;

	if H[1] <> H[Length(H)] then
		# the first vertex is not the same as the last vertex
		return false;
	elif Length(H) <> OrderGraph(X) + 1 then
		# the cycle does not have the correct length
		return false;
	elif not(ForAll(
                [1 .. Length(H)-1],
                x -> ( (H{[x,x+1]} in EdgesOfX) or (H{[x,x+1]} in AdditionalEdges) )
                ))
    then
        # there are consecutive vertices that are not joined by an edge
        return false;
    elif AsSSortedList(H) <> AsSSortedList(VertexNames(X)) then
		# the cycle does not pass through every vertex
		return false;
	elif not(ForAll( RequiredEdges, edge -> 
				(edge in List([1..Length(H)-1], i -> H{[i,i+1]}))
				or ( (Reversed(edge) in RequiredEdges) and (Reversed(edge) in List([1..Length(H)-1], i -> H{[i,i+1]})))))
        then
		# the hamiltonian cycle does not include all of the required edges
		return false;
	else
		# H passed all the tests, so it is a valid hamiltonian cycle
		return true;
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

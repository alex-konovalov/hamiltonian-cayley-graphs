# AbelianCharacters.gap

AbelianCharacters := function(G)
    # This function returns a list of all the 1-dimensional characters of the group G.
    
    local A, pA, CharacterList, chi, ConvertMappingToFunction;
    
    # Let A be the abelianization of G and 
    # let pA be the natural homomorphism from G to A
    A := FactorGroup(G, CommutatorSubgroup(G,G) );
    pA := NaturalHomomorphism(A);
    
    # 1-dimensional characters of G are the same as irreducible representations of
    # the abelianization A of G: each irreducible representation chi of A yields a 
    # 1-dimensional character of G after composing it with the homomorphism pA from 
    # G to A.
    CharacterList := List(IrreducibleRepresentations(A), chi -> CompositionMapping(chi,pA) );

    # CharacterList is now a list of "mappings", but they need to be converted to functions.
    # Also, the value of the character is a (one-by-one) matrix representation, but we
    # want a scalar, so we strip off the matrix brackets by returning the (1,1) entry of
    # the matrix.
    ConvertMappingToFunction := function(aMatrixMap)
        return x -> ImageElm(aMatrixMap,x)[1][1];
    end;
    return List(CharacterList, ConvertMappingToFunction);
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




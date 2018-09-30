# CheckForIrredundantSpecialCase.gap

CheckForIrredundantSpecialCase := function(GBar, SBar, twist, k, GapId)
	# check whether we are in a special case that is covered by Lemma 3.1
	# (or some other result), and therefore does not need to be considered.  
	# If so, return the reason (as a string). Otherwise, return fail.	``

	if IsAbelian(GBar) and ForAll(SBar, s -> twist(s)=1) then
		return "G is abelian";

	elif (Length(SBar) = 1) and (twist(SBar[1]) <> 1)
	then return "S cannot generate G because nonabelian groups are not cyclic";

	elif
		Length(SBar) = 2 
			and 
		(
			(Order(SBar[1]) = 2 and Order(SBar[2]) = 3 and twist(SBar[2]) = 1) 
				or 
			(Order(SBar[2]) = 2 and Order(SBar[1]) = 3 and twist(SBar[1]) = 1)
		)
	then return "Lemma 3.1(1) applies";

	elif
		List(SBar, Order) = [3,3]
			and
		List(SBar, twist) = [1,1]
			and
		[k, GapId] = IdSmallGroup(AlternatingGroup(4))
	then return "Lemma 3.1(2) applies";

	elif 
		Length(SBar) = 2
			and
		AsSortedList( List(SBar, Order)) = [2, k/4] 
			and
		IsOddInt(k/4)
			and
		List(SBar, twist) = [1,1]
			and
		IdSmallGroup( CommutatorSubgroup(GBar,GBar) ) = IdSmallGroup(DihedralGroup(4))
	then return "Lemma 3.1(3) applies";

	elif ForAny( SBar, s -> (s in Centre(GBar) and Order(s) = 2 and twist(s) = 1) )
	then return "Lemma 3.1(4) applies";
									
	elif ForAny( SBar, s -> (IsPrimeInt(Order(s^2)) and s in Centre(GBar) and twist(s) = -1) )
	then return "Lemma 3.1(5) applies";
	
	elif 
		Length(SBar) = 2
			and
		IsEvenInt(k)
			and
		(
			(Order(SBar[1]) = 2 and twist(SBar[1]) = -1 
				and Order(SBar[2]) = k/2 and twist(SBar[2]) = 1)
			or
			(Order(SBar[2]) = 2 and twist(SBar[2]) = -1 
				and Order(SBar[1]) = k/2 and twist(SBar[1]) = 1)
		)
			and
		[k,GapId] = IdSmallGroup(DihedralGroup(k))
	then return "Lemma 3.1(6) applies";

	elif
		Order(GBar) = 20
			and
		Order(CommutatorSubgroup(GBar,GBar)) = 5
			and
		AsSortedList(List(SBar, Order)) = [2,4]
			and
		List(SBar, twist) = [1,1]
			and
		CommutatorSubgroup(GBar,GBar).1 * SBar[1] <> SBar[1] * CommutatorSubgroup(GBar,GBar).1
			and
		CommutatorSubgroup(GBar,GBar).1 * SBar[2] <> SBar[2] * CommutatorSubgroup(GBar,GBar).1
	then return "Lemma 3.1(7) applies (with q = 5)";
	else 
		return fail;
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

# FindNonzeroDet.gap

Read("LKH.gap");
Read("CyclotomicToFp.gap");
Read("Det.gap");

FindNonzeroDetRoutine := function(AllVoltages, MatSize, MatSoFar, RowsUsed, WhereToStart)
    local LastRowUsed, RowToTry, NextMatrix;
    # recursively try to find MatSize linearly independent rows, by greedily choosing
    # the first row that is linearly independent from all previous rows
    
    LastRowUsed := Maximum(Concatenation(RowsUsed, [WhereToStart-1]));
    if Length(MatSoFar) = MatSize then
        if Det(MatSoFar) = 0 then
            return fail;
        else
            return [MatSoFar, RowsUsed];
        fi;
    else
        if LastRowUsed < Length(AllVoltages) then
            for RowToTry in [(LastRowUsed + 1)..Length(AllVoltages)] do
                NextMatrix := Concatenation(MatSoFar, [AllVoltages[RowToTry]]);
                if RankMat(NextMatrix) = Length(NextMatrix) then
                    return FindNonzeroDetRoutine(AllVoltages, MatSize, NextMatrix, Concatenation(RowsUsed, [RowToTry]), WhereToStart);
                fi;
            od;
        fi;
        return fail;
    fi;
end;

FindNonzeroDet := function(AllVoltages)
	# call FindNonzeroDetRoutine to find rows that provide a nonzero determinant.
	# If we are working in the cyclotomic numbers (ie., not in a finite field)
	# and FindNonzeroDetRoutine provides a determinant whose norm has a huge prime divisor, 
	# then we tell it to try again, leaving out the first row that used last time.
	
	local WhereToStart, MaxPrimeAllowed, MatSize, MatAndRows;
    MaxPrimeAllowed := 10^5;
    
	WhereToStart := 1;
	MatSize := Length(AllVoltages[1]);
	
	while WhereToStart + MatSize - 1 <= Length(AllVoltages) do
		MatAndRows := FindNonzeroDetRoutine(AllVoltages, MatSize, [], [], WhereToStart);
		if (MatAndRows = fail) then
			return fail;
		elif Length(MatAndRows[1]) = 0 or IsFFE(MatAndRows[1][1][1]) or Maximum(Concatenation([1], PrimeDivisors(Norm(Det(MatAndRows[1]))))) < MaxPrimeAllowed then
			return MatAndRows[2];
		else
			WhereToStart := MatAndRows[2][1] + 1 ; # start looking again, past the last one
		fi;
	od;
	return fail;
end;

FindNonzeroDetModp := function(AllVoltages, p)
	# instead of taking a norm, we reduce all voltages modulo p, and look for a 
	# nonzero determinant
	local AllVoltagesModp;
	AllVoltagesModp := List(AllVoltages, VoltageList -> List(VoltageList, voltage -> CyclotomicToFp(voltage,p)));
	return FindNonzeroDet(AllVoltagesModp);
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

# The function HamiltonianCycles(X, RequiredEdges) returns a list of all the hamiltonian 
# cycles in the graph X that contain all of the edges in the list RequiredEdges.  (X must 
# be a GRAPE graph.)  It is intended only for fairly small graphs, with small valence and
# no more than a few dozen vertices.
#
# Actually, the function returns an iterator, rather than a list, because there may
# be exponentially many hamiltonian cycles.
#
# Each hamiltonian cycle is returned as a list of the names of the vertices in X,
# in the order that they are visited.

LoadPackage("grape");

RemoveEdge := function(X, u, v)
	# remove the edge uv from consideration
	Remove(X.adjacencies[u], Position(X.adjacencies[u], v));
	Remove(X.adjacencies[v], Position(X.adjacencies[v], u));
end;

UseEdge := function(X, u, v)
	# add edge uv to the hamiltonian cycle we are building, and return true
	# return false if this results in a premature cycle or a cycle vertex of valence > 2
	local pathEndu, pathEndv;
	if 
		( (X.pathEnd[u] <> v) or (X.pathLength[u] = X.order - 1) ) 
			and 
		(Length(X.cycleEdges[u]) < 2) and (Length(X.cycleEdges[v]) < 2) 
	then
		Add( X.cycleEdges[u], v );
		Add( X.cycleEdges[v], u );
		pathEndu := X.pathEnd[u];
		pathEndv := X.pathEnd[v];
		X.pathEnd[pathEndu] := pathEndv;
		X.pathEnd[pathEndv] := pathEndu;
		X.pathLength[pathEndu] := X.pathLength[u] + X.pathLength[v] + 1;
		X.pathLength[pathEndv] := X.pathLength[pathEndu];
		return true;
	else
		return false;
	fi;
end;

IsConnected := function(X)
	# this is a replacement for grape's built-in IsConnectedGraph function,
	# because the built-in function is too slow on our graphs
	local reachOut, canReach;
	
	reachOut := function(X, u, canReach)
		# recursively set canReach to true for all vertices that can be reached from u
		local v;
		canReach[u] := true;
		for v in X.adjacencies[u] do
			if canReach[v] = false then
				reachOut(X, v, canReach);
			fi;
		od;
	end;

	canReach := List([1..X.order], x -> false);
	reachOut(X,1,canReach);
	return ForAll(canReach, x -> x);
end;

FindHamCyc := function(X, Data)
	local n, V, W, u, v, i, j, vert, hasChanged, valid, connected, canReach, 
		LengthOfCycle, previousVertex,
		NextGraph, result, valenceOfu, cycleValenceOfu;
	
	n := X.order;
	
	# add/remove all forced edges
	hasChanged := true;
	valid := true;
	while hasChanged and valid do
		hasChanged := false;
		for u in [1..n] do
			if valid then
				valenceOfu := Length(X.adjacencies[u]);
				cycleValenceOfu := Length(X.cycleEdges[u]);
				if (valenceOfu < 2) then
					valid := false;
					break;
				elif (valenceOfu = 2) and (cycleValenceOfu < 2) then
					hasChanged := true;
					for v in X.adjacencies[u] do
						if not (v in X.cycleEdges[u]) then
							if not(UseEdge(X, u, v)) then
								valid := false;
							fi;
						fi;
					od;
				elif (cycleValenceOfu = 2) and (valenceOfu > 2) then
					hasChanged := true;
					for v in X.adjacencies[u] do
						if not (v in X.cycleEdges[u]) then
							RemoveEdge(X, u, v);
						fi;
					od;
				elif (cycleValenceOfu = 1) and (X.pathEnd[u] in X.adjacencies[u]) and not(X.pathEnd[u] in X.cycleEdges[u]) and (X.pathLength[u] <> X.order - 1) then
					hasChanged := true;
					RemoveEdge(X, u, X.pathEnd[u]);
				fi;
			fi;
		od;
	od;

	# check that the graph is still connected
	if valid then
		connected := IsConnected(X);
		if not(connected) then
			valid := false;
		fi;
	fi;
	
	v := fail;
	if valid then
		# find an edge to try
		u := 0;
		while (v = fail) and (u < n) do
			u := u + 1;
			v := First(X.adjacencies[u], x -> not(x in X.cycleEdges[u]));
		od;
	fi;
	
	if v <> fail then
	
		# we will try using the edge, but first make a copy of the graph without this edge
		NextGraph := StructuralCopy(X);
		RemoveEdge(NextGraph, u, v);
		
		if UseEdge(X, u,v) then
			result := FindHamCyc(X, Data);
			if result <> fail then
				Add(Data, NextGraph, 1);
				return result;
			else
				return FindHamCyc(NextGraph, Data);
			fi;
		else
			return FindHamCyc(NextGraph, Data);
		fi;
	else
		# if the "cycle" is valid and connected and every vertex has valence 2 then
		# we have a hamiltonian cycle !!!
		if valid and connected and ForAll(X.cycleEdges, x -> Length(x) = 2) then 
			return X;	
		else
			return fail;
		fi;
	fi;
end;




NoMoreHamCycs := function(iter)
	local X, Data, result;

	Data := StructuralCopy(iter!.stack);
	if not(iter!.UsedThisHC) then
		Add(Data, iter!.graph, 1);
	fi;
	iter!.UsedThisHC := false;
	
	while Data <> [] do
		X := Data[1];
		Remove(Data,1);
		result := FindHamCyc(X, Data);
		if result <> fail then
			iter!.graph := result;
			iter!.stack := StructuralCopy(Data);
			return false;
		fi;
	od;
	return true;
end;
	

HamiltonianCycles := function (X, RequiredEdges)
	local Xcopy, v, e;
	
	# make a copy of X that can be modified by the algorithm
	# (but without multiple edges or loops)	
	Xcopy := rec(
		order := X.order,
		adjacencies := List([1..X.order], 
			function(v) 
				return 
					Filtered(
						# AsSSortedList removes repeated edges
						AsSSortedList(Adjacency(X,v)), 
						# filtering by x -> x <> v removes loops
						x -> x <> v
					); 
			end
		),
		cycleEdges := List([1..X.order], x -> []),
		pathEnd := [1..X.order],
		pathLength := List([1..X.order], x -> 0)
		);
	# put all required edges into the hamiltonian cycle
	for e in RequiredEdges do
		UseEdge(Xcopy, Position(VertexNames(X),e[1]), Position(VertexNames(X),e[2]) );
	od;	
	
	return IteratorByFunctions( rec(
		
		NextIterator := function(iter)
			local i, u, v, Hedges, H;
			if not(NoMoreHamCycs(iter)) then
				Hedges := iter!.graph.cycleEdges;
				H := [ VertexName(iter!.grape, 1) ];
				u := 1;
				v := Hedges[1][1];
				for i in [1..iter!.grape.order] do
					Add(H, VertexName(iter!.grape, v));
					if Hedges[v][1] = u then
						u := v;
						v := Hedges[v][2];
					else
						u := v;
						v := Hedges[v][1];
					fi;
				od;
			
				iter!.UsedThisHC := true;
				return H;
			fi;
		end,

		IsDoneIterator := NoMoreHamCycs,

		ShallowCopy := function(iter)
				return rec(
					NextIterator := iter!.NextIterator,
					IsDoneIterator := iter!.IsDoneIterator,
					ShallowCopy := iter!.ShallowCopy,
					PrintObj := iter!.PrintObj,
					grape := ShallowCopy(iter!.grape),
					graph := ShallowCopy(iter!.graph),
					stack := ShallowCopy(iter!.stack),
					UsedThisHC := iter!.UsedThisHC
					);
				end,

		PrintObj := function(iter)
			Print(
				"Iterator for hamiltonian cycles in a graph with ", 
				iter!.grape.order, " vertices."
			);
		end,

		grape := X, # keep a copy of the original graph X
		
		graph := Xcopy,
		
		stack := [],
		
		UsedThisHC := false
	));
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

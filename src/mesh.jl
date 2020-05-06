function order_triangles!(t, tnbor)
	# do breadth first search to make normals consistent

	visited = falses(size(t,2))
	q = Array{Int64}(undef,0)
	S = Set{Tuple{Int64,Int64}}()

	push!(q, 1)

	while !isempty(q)
		i = popfirst!(q)

		if visited[i]
			continue
		end

		visited[i] = true

		i1 = t[1,i]; i2 = t[2,i]; i3 = t[3,i];

		if (i1,i2) ∈ S || (i2,i3) ∈ S || (i3,i1) ∈ S
			# need to reverse chirality of triangle
			swap = i2
			i2 = i3
			i3 = swap
			t[2,i] = i2; t[3,i] = i3; # update triangle indices
		end
		
		# add edges to set
		push!(S, (i1,i2))
		push!(S, (i2,i3))
		push!(S, (i3,i1))

		# add triangle neighbors to queue
		if (j = tnbor[1,i]) > 0; push!(q, j); end
		if (j = tnbor[2,i]) > 0; push!(q, j); end
		if (j = tnbor[3,i]) > 0; push!(q, j); end
	end
end

function calculate_neighbors(t, Np)
	# given a list of indices of triangles, t, returns a list of neighbors for each particle
	
	maxN = 10

	nbor_count = zeros(Int64, Np) # length of filled neighbors
	nbor = zeros(Int64, maxN, Np) # assumes no more than 10 neighbors
	tmem = zeros(Int64, maxN, Np) # triangles that each point is a member of
	tmem_count = zeros(Int64, Np)

	function check_count(nind, tind)
		if (nind > 0 && nbor_count[nind] > maxN) || (tind > 0 && tmem_count[tind] > maxN)
			new_nbor = zeros(Int64, 2*maxN, Np)
			new_tmem = zeros(Int64, 2*maxN, Np)
			new_nbor[1:maxN,:] .= nbor
			new_tmem[1:maxN,:] .= tmem
			nbor = new_nbor
			tmem = new_tmem
			maxN *= 2
		end
	end
	
	for i = 1:size(t,2)
		i1 = t[1,i]; i2 = t[2,i]; i3 = t[3,i];
		
		# iterate through neighbors
		if !any(nbor[:, i1] .== i2)
			nbor_count[i1] += 1
			check_count(i1, 0)
			nbor[nbor_count[i1], i1] = i2
		end
		if !any(nbor[:, i1] .== i3)
			nbor_count[i1] += 1
			check_count(i1, 0)
			nbor[nbor_count[i1], i1] = i3
		end
		if !any(nbor[:, i2] .== i1)
			nbor_count[i2] += 1
			check_count(i2, 0)
			nbor[nbor_count[i2], i2] = i1
		end
		if !any(nbor[:, i2] .== i3)
			nbor_count[i2] += 1
			check_count(i2, 0)
			nbor[nbor_count[i2], i2] = i3
		end
		if !any(nbor[:, i3] .== i1)
			nbor_count[i3] += 1
			check_count(i3, 0)
			nbor[nbor_count[i3], i3] = i1
		end
		if !any(nbor[:, i3] .== i2)
			nbor_count[i3] += 1
			check_count(i3, 0)
			nbor[nbor_count[i3], i3] = i2
		end
		
		# add triangle to triangle members list
		if !any(tmem[:, i1] .== i)
			tmem_count[i1] += 1
			check_count(0, i1)
			tmem[tmem_count[i1], i1] = i
		end
		if !any(tmem[:, i2] .== i)
			tmem_count[i2] += 1
			check_count(0, i2)
			tmem[tmem_count[i2], i2] = i
		end
		if !any(tmem[:, i3] .== i)
			tmem_count[i3] += 1
			check_count(0, i3)
			tmem[tmem_count[i3], i3] = i
		end
	end
	
	return nbor, nbor_count, tmem, tmem_count
end

function triangle_connectivity(t)
	# calculates triangle neighbors given a list of indices of triangle points

	Nt = size(t,2)
	tnbor_count = zeros(Int64, Nt)
	tnbor = zeros(Int64, 3, Nt)

	D = Dict{Tuple{Int64,Int64}, Int64}()

	for t1 = 1:Nt
		for (j1,j2) in [(1,2), (2,3), (3,1)]
			i1 = t[j1,t1]; i2 = t[j2,t1];

			edge = (min(i1,i2), max(i1,i2))
			if haskey(D, edge)
				# add neighbors since edge was shared
				t2 = D[edge]
				tnbor_count[t1] += 1
				tnbor_count[t2] += 1
				tnbor[tnbor_count[t1], t1] = t2
				tnbor[tnbor_count[t2], t2] = t1
			else
				# store triangle in dictionary
				push!(D, edge => t1)
			end
		end
	end

	return tnbor
end

function find_boundary_vertices(nbor_count, tmem_count)
	#= Finds the vertices that lie on the boundary of the mesh in O(Np) time
	A vertex is on the boundary if the number of neighbors != number of triangles
	
	INPUT
	nbor_count - number of neighbors for each vertex
	tmem_count - number of triangles to which a vertex belongs
	
	
	OUTPUT
	isbndry - a boolean vector
	=#
	
	isbndry = (nbor_count .!= tmem_count)
	
	return isbndry
end

function find_boundary_edges(isbndry, t)
	#= find the path of vertices that lie on the boundary of a surface
	
	INPUT
	isbndry - vector of booleans for border vertices
	t
	
	OUTPUT
	bndry_path - vector containing indices of edges on the boundary of a surface
		separate boundaries are separated by 0 entries
		one boundary is always terminated with a 0 entry
	=#
	
	Np = length(isbndry)
	
	# set of edges that have only been visited once (i.e., on boundary)
	edges = Set{Tuple{Int64,Int64}}()
	
	# loop through all edges of triangles
	for i = 1:size(t,2)
		for (j1,j2) in [(1,2), (2,3), (3,1)]
			i1 = t[j1,i]; i2 = t[j2,i];

			edge = (min(i1,i2), max(i1,i2))
			if edge ∈ edges
				delete!(edges, edge)
			else
				push!(edges, edge)
			end

		end
	end

	# dictionary to map points to their neighbors on the boundary
	map_to_neigh = Dict{Int64, Tuple{Int64,Int64}}()

	for (i1,i2) in edges
		if haskey(map_to_neigh, i1)
			push!(map_to_neigh, i1 => (map_to_neigh[i1][1], i2)) # i2 is the second neighbor found
		else
			push!(map_to_neigh, i1 => (i2, 0))
		end

		# do the same for i2
		if haskey(map_to_neigh, i2)
			push!(map_to_neigh, i2 => (map_to_neigh[i2][1], i1)) # i2 is the second neighbor found
		else
			push!(map_to_neigh, i2 => (i1, 0))
		end
	end
	
	# traverse all of the boundary edges
	bndry_path = zeros(Int64, 0)

	# get neighbor of ind and remove ind from neighbor's neighbors
	function next_ind(ind)
		(i1,i2) = map_to_neigh[ind]

		if i2 == 0
			delete!(map_to_neigh, ind)
		else
			map_to_neigh[ind] = (i2, 0)
		end

		# remove recripocal neighbor mapping from i1
		(i3,i4) = map_to_neigh[i1]
		@assert (i3 == ind || i4 == ind)

		if i4 == 0
			delete!(map_to_neigh, i1)
		elseif i4 == ind
			map_to_neigh[i1] = (i3, 0)
		else
			map_to_neigh[i1] = (i4, 0)
		end

		return i1
	end
	
	# iterate through map of neighbors until empty
	while !isempty(map_to_neigh)
		ind = first(keys(map_to_neigh))
		push!(bndry_path, ind)

		# keep traversing neighbors until loop ends
		while haskey(map_to_neigh, ind)
			ind = next_ind(ind)
			push!(bndry_path, ind)
		end

		# append loop with a 0
		push!(bndry_path, 0)
	end
	   
	return bndry_path
end

function contour_length(p, bndry_path)
	#= calculates the countour length of the boundaries
	
	INPUT
	p - a stacked vector of coordinates [x1 y1 z1 x2 y2 z2 ...]
	bndry_path - edges of boundaries
	
	OUTPUT
	val - total contour length
	=#
   
	val = zero(eltype(p))
	
	for i = 1:length(bndry_path)-1
		if bndry_path[i] == 0 || bndry_path[i+1] == 0
			continue
		end
		
		i1 = bndry_path[i]; i2 = bndry_path[i+1];
		
		# calculate bond vectors
		v12 = p[3*i2-2:3*i2] .- p[3*i1-2:3*i1]
		
		val += norm(v12)
	end
	
	return val
end

function area_volume(p, t)
	#= calculates the surface area of a mesh
	calculates the volume enclosed by a mesh using signed volumes
	t needs to have a consistent ordering (so normals all face the same way)
	
	INPUT
	p - a stacked vector of coordinates [x1 y1 z1 x2 y2 z2 ...]
	t - indices of triangles
	
	OUTPUT
	A - area
	V - volume
	=#
	
	A = zero(eltype(p))
	V = zero(eltype(p))
	
	for i = 1:size(t,2)
		i1 = t[1,i]; i2 = t[2,i]; i3 = t[3,i];
		
		# calculate bond vectors
		v1 = p[3*i1-2:3*i1]
		v2 = p[3*i2-2:3*i2]
		v3 = p[3*i3-2:3*i3]
		
		# determine the side of the triangle on which the origin lies
		sgn = 1;
		n = cross(v1 .- v2, v3 .- v2)
		if dot(-v1, n) < 0
			sgn = -1
		end
		
		A += norm(n)/2
		
		# divide by 6 since it's the volume of a tetrahedron
		V += sgn * (v1[1]*(v2[2]*v3[3] - v2[3]*v3[2]) - v1[2]*(v2[1]*v3[3] - v2[3]*v3[1]) + v1[3]*(v2[1]*v3[2] - v2[2]*v3[1]))/6
	end
	
	# abs because it could be the wrong sign
	return A, abs(V)
end

function area(p, t)
	#= calculates the surface area of a mesh
	
	INPUT
	p - a stacked vector of coordinates [x1 y1 z1 x2 y2 z2 ...]
	t - indices of triangles
	
	OUTPUT
	A - area
	=#
	
	A = zero(eltype(p))
	
	for i = 1:size(t,2)
		i1 = t[1,i]; i2 = t[2,i]; i3 = t[3,i];
		
		# calculate bond vectors
		v1 = p[3*i1-2:3*i1]
		v2 = p[3*i2-2:3*i2]
		v3 = p[3*i3-2:3*i3]
		
		A += norm(cross(v1 .- v2, v3 .- v2))/2
	end
	
	# abs because it could be the wrong sign
	return A
end

function areas(p::Vector, t)
	#= calculates the surface area of each triangle in a mesh
	
	INPUT
	p - a stacked vector of coordinates [x1 y1 z1 x2 y2 z2 ...]
	t - indices of triangles
	
	OUTPUT
	A - area
	=#
	
	A = zeros(size(t,2))
	
	for i = 1:size(t,2)
		i1 = t[1,i]; i2 = t[2,i]; i3 = t[3,i];
		
		# calculate bond vectors
		v1 = p[3*i1-2:3*i1]
		v2 = p[3*i2-2:3*i2]
		v3 = p[3*i3-2:3*i3]
		
		A[i] = norm(cross(v1 .- v2, v3 .- v2))/2
	end
	
	# abs because it could be the wrong sign
	return A
end

function normals(p::Vector, t)
	# returns the normals of all of the triangles
	N = zeros(3, size(t,2))

	for t1 = 1:size(t,2)
		i1 = t[1,t1]; i2 = t[2,t1]; i3 = t[3,t1];

		# calculate bond vectors and normals for each triangle
		v12 = p[3*i2-2:3*i2] .- p[3*i1-2:3*i1]
		v13 = p[3*i3-2:3*i3] .- p[3*i1-2:3*i1]

		n = cross(v12, v13)
		normalize!(n)

		N[:,t1] .= n
	end

	return N
end




# ============================================ REMESHING =================================================






function subdivide(x::Vector, mi::MeshInfo)
	#= subdivides a mesh by splitting each triangle into 4 triangles (like the icosahedral algorithm)

	INPUT
	x - vertices of mesh
	mi - MeshInfo object

	OUTPUT
	xnew, minew - new vertex list and MeshInfo object
	=#

	Np = length(x) ÷ 3

	# new points and triangle indices
	xnew = copy(x)
	resize!(xnew, length(x) + 3*(sum(mi.nbor_count) ÷ 2))
	tnew = similar(mi.t, 3, 4*size(mi.t, 2))

	t = mi.t

	# set of edges
	D = Dict{Tuple{Int64,Int64}, Int64}()
	c = Np + 1

	# loop through edges to add new points
	for i = 1:Np
		for j = 1:mi.nbor_count[i]
			i1 = i; i2 = mi.nbor[j,i]

			# add point to edge if it hasn't been added already
			if !haskey(D, (min(i1,i2), max(i1,i2)))
				xnew[3*c-2:3*c] .= (x[3*i1-2:3*i1] .+ x[3*i2-2:3*i2])/2

				# add edge to dictionary along with new point
				push!(D, (min(i1,i2), max(i1,i2)) => c)
				c += 1
			end
		end
	end

	# construct new triangles
	for i = 1:size(t, 2)
		i1 = t[1,i]; i2 = t[2,i]; i3 = t[3,i];

		tnew[:, 4*i-3] .= [i1, D[(min(i1,i2), max(i1,i2))], D[(min(i1,i3), max(i1,i3))]]
		tnew[:, 4*i-2] .= [i2, D[(min(i1,i2), max(i1,i2))], D[(min(i2,i3), max(i2,i3))]]
		tnew[:, 4*i-1] .= [i3, D[(min(i1,i3), max(i1,i3))], D[(min(i2,i3), max(i2,i3))]]
		tnew[:, 4*i-0] .= [D[(min(i1,i2), max(i1,i2))], D[(min(i1,i3), max(i1,i3))], D[(min(i2,i3), max(i2,i3))]]
	end

	mi_new = MeshInfo(tnew)

	# ensure triangles have consistent ordering
	order_triangles!(tnew, mi_new.tnbor)

	return xnew, mi_new
end


function mean_edge_length(x::Vector, mi::MeshInfo)
	# calculates the mean edge length of a mesh
	l = 0.0
	Nedges_x2 = sum(mi.nbor_count)

	for i = 1:size(mi.nbor,2)
		for j = 1:mi.nbor_count[i]
			i1 = i; i2 = mi.nbor[j,i]

			l += norm(x[3*i2-2:3*i2] .- x[3*i1-2:3*i1])
		end
	end

	return l / Nedges_x2
end

function construct_edge_map(t::Array{Int64,2})
	#= constructs a map from edges (stored as tuples of ordered vertex indices) to
	triangles that they touch (also tuples, where the second integer is 0 if an edge
	only touches one triangle)

	INPUT
	t - triangle indices

	OUTPUT
	edge_map - dictionary object that maps tuples to tuples
	=#

	edge_map = Dict{Tuple{Int64,Int64}, Tuple{Int64,Int64}}()

	@inbounds for t1 = 1:size(t,2)
		i1 = t[1,t1]; i2 = t[2,t1]; i3 = t[3,t1];

		edge1 = (min(i1,i2), max(i1,i2))
		edge2 = (min(i2,i3), max(i2,i3))
		edge3 = (min(i3,i1), max(i3,i1))

		for edge in [edge1, edge2, edge3]
			if haskey(edge_map, edge)
				push!(edge_map, edge => (edge_map[edge][1], t1))
			else
				push!(edge_map, edge => (t1, 0))
			end
		end
	end

	return edge_map
end

function collapse_short_edges(x::Vector, mi::MeshInfo, l::Float64)
	collapse_short_edges(x, mi, p -> l)
end


function collapse_short_edges(x::Vector, mi::MeshInfo, l::Function)
	#= Collapse edges shorter than a factor less than the mean length in a greedy manner

	INPUT
	x - vector of positions
	mi - the mesh
	l - collapses edges of length less than l

	OUTPUT
	new_x - new vector of positions
	new_mi - new MeshInfo object
	=#

	xc = copy(x)
	Np = length(x) ÷ 3

	# construct neighbor data
	# (using an array of arrays for vertex neighbors for easier unions)
	nbor = [mi.nbor[1:mi.nbor_count[i],i] for i in 1:Np]
	tnbor = copy(mi.tnbor)

	edge_map = construct_edge_map(mi.t)

	# maps deleted points to their equivalence classes
	equiv_map = Dict(i => i for i in 1:Np)

	# map of old points to the new collapsed points (when collapsing an edge)
	equiv_collapse = [[i] for i in 1:Np]

	# put all edges into queue (adding them only once)
	q = collect(keys(edge_map))

	# return triangles neighboring t1 that aren't t2
	@inline function get_other_tris(t1, t2)
		t11, t12, t13 = tnbor[1:3, t1]
		t11 == t2 && return t12, t13
		t12 == t2 && return t11, t13
		return t11, t12
	end

	# replace neighboring triangle t2 of t1 with t3
	@inline function replace_tri(t1, t2, t3)
		if t1 == 0
			return
		end

		@assert t2 != 0
		@assert t2 ∈ tnbor[:,t1]

		tnbor[1,t1] == t2 && (tnbor[1,t1] = t3)
		tnbor[2,t1] == t2 && (tnbor[2,t1] = t3)
		tnbor[3,t1] == t2 && (tnbor[3,t1] = t3)
	end

	# get other point in triangle tind that's not in edge
	# needs to use equiv_map since mi.t contains old, non merged indices!
	@inline function get_other_point_in_tri(tind, edge)
		eq1 = equiv_map[mi.t[1,tind]]
		eq2 = equiv_map[mi.t[2,tind]]
		eq3 = equiv_map[mi.t[3,tind]]

		if eq1 != edge[1] && eq1 != edge[2]
			return eq1
		elseif eq2 != edge[1] && eq2 != edge[2]
			return eq2
		else
			return eq3
		end
	end

	i2e = (a,b) -> (min(a,b), max(a,b))

	# find edges to delete
	while !isempty(q)
		edge = popfirst!(q)

		# skip if edge doesn't exist anymore
		if !haskey(edge_map, edge)
			continue
		end

		(i1, i2) = edge
		@assert i1 < i2

		# don't collapse edges if they're touching a vertex that has already been collapsed
		if length(equiv_collapse[i1]) > 1 || length(equiv_collapse[i2]) > 1
			continue
		end

		# skip edge if touching boundary
		if mi.isbndry[i1] || mi.isbndry[i2]
			continue
		end

		@assert equiv_map[i1] == i1
		@assert equiv_map[i2] == i2
		@assert i1 != i2

		p1 = (xc[3*i1-2], xc[3*i1-1], xc[3*i1])
		p2 = (xc[3*i2-2], xc[3*i2-1], xc[3*i2])
		mean_point = (p1 .+ p2)./2

		edge_length = norm(p2 .- p1)

		# don't collapse edges on the boundary
		if edge_length < l(mean_point)

			# mat"figure(3); hold on; plot3([$x(3*$i2-2) $x(3*$i1-2)], [$x(3*$i2-1) $x(3*$i1-1)], [$x(3*$i2) $x(3*$i1)], 'r-', 'linewidth', 2); hold off"

			@assert i1 != i2

			# update triangle neighbors
			t1, t2 = edge_map[edge]
			@assert t1 != 0 && t2 != 0	# should both be nonzero since boundary is prohibited

			t11, t12 = get_other_tris(t1, t2)
			t21, t22 = get_other_tris(t2, t1)

			# ensure that t11 and t12 and t21 and t22 aren't neighbors themselves (otherwise would create degeneracy)
			if t12 == tnbor[1,t11] || t12 == tnbor[2,t11] || t12 == tnbor[3,t11]
				continue
			elseif t22 == tnbor[1,t21] || t22 == tnbor[2,t21] || t22 == tnbor[3,t21]
				continue
			end

			# ensure that the neighboring triangles don't equal each other
			if t11 == t21 || t12 == t22 || t11 == t12 || t21 == t22 || t21 == t12 || t11 == t22
				continue
			end

			replace_tri(t11, t1, t12)
			replace_tri(t12, t1, t11)
			replace_tri(t21, t2, t22)
			replace_tri(t22, t2, t21)
			tnbor[:, t1] .= 0		# t1 and t2 no longer exist!
			tnbor[:, t2] .= 0

			# update edge map (which should only contain equiv class references based on update procedure)
			t1_other = get_other_point_in_tri(t1, edge)
			t2_other = get_other_point_in_tri(t2, edge)

			for i3 in nbor[i2]									# delete occurences of i2 in edge_map
				ts_touching = pop!(edge_map, i2e(i2, i3))

				if i3 != i1 && i3 != t1_other && i3 != t2_other
					push!(edge_map, i2e(i1, i3) => ts_touching) # update edge with i1 instead of i2
					push!(q, i2e(i1, i3))						# push updated edge to queue
				end

				deleteat!(nbor[i3], findfirst(isequal(i2), nbor[i3]))	# update i2's neighbors
				i3 != i1 && union!(nbor[i3], i1)
			end

			edge_map[i2e(i1, t1_other)] = (max(t11, t12), min(t11, t12))	# update triangles in other 2 edges
			edge_map[i2e(i1, t2_other)] = (max(t21, t22), min(t21, t22))
			
			# update point neighbors and push new edges
			# although these operations require looping through all neighbors, this should be more efficient than a set
			deleteat!(nbor[i2], findfirst(isequal(i1), nbor[i2]))	# remove occurence of i1
			union!(nbor[i1], nbor[i2])
			resize!(nbor[i2], 0)

			# merge classes and udpate equiv_map
			equiv_class1 = equiv_collapse[i1]
			equiv_class2 = equiv_collapse[i2]

			for el in equiv_class2
				push!(equiv_map, el => i1)
			end

			append!(equiv_class1, equiv_class2)		# move all elements in 2 to 1
			resize!(equiv_class2, 0)

			# finally update point
			xc[3*i1-2:3*i1] .= mean_point

		end
	end

	# ------------- construct new set of points and new triangulation ---------------
	new_t = zeros(Int64, 0)
	new_x = zeros(0)

	# to do remapping of indices from old points to new points
	new_point_map = Dict{Int64, Int64}()

	# returns the new index of the point
	# if it has not been used yet, it creates the new mapping and adds the point to new_x
	function get_new_index(ind)
		# append point to new_x if it hasn't been mapped yet
		return get(new_point_map, ind) do
			append!(new_x, xc[3*ind-2:3*ind])
			newind = length(new_x) ÷ 3
			push!(new_point_map, ind => newind)
			newind
		end
	end

	# loop through all triangles and reconstruct new triangles and points
	for t1 = 1:size(mi.t,2)
		# skip if triangle has been deleted
		if tnbor[1,t1] == 0
			continue
		end

		i1 = equiv_map[mi.t[1,t1]]; i2 = equiv_map[mi.t[2,t1]]; i3 = equiv_map[mi.t[3,t1]];

		new_i1 = get_new_index(i1);
		new_i2 = get_new_index(i2);
		new_i3 = get_new_index(i3);

		push!(new_t, new_i1, new_i2, new_i3)
	end

	t_as_mat = reshape(new_t, (3, length(new_t) ÷ 3))

	# new_x = copy(x) # ====================================
	# add_text_labels(new_x)
	# @show t_as_mat'

	# new_tnbor = triangle_connectivity(t_as_mat)
	# order_triangles!(t_as_mat, new_tnbor)
	try
		new_tnbor = triangle_connectivity(t_as_mat)
		order_triangles!(t_as_mat, new_tnbor)
	catch err
		# trisurf(new_x, t_as_mat, 2)
		# # @show t_as_mat'
		# @show err.i[2]
		# mat"figure(2); hold on;"
		# trisurf(new_x, t_as_mat[:,err.i[2]], 2, edgecolor="r", facecolor="none")
		# mat"figure(2); hold off;"
		# throw(err)
		@assert false "Mesh degenerated"
	end

	return new_x, MeshInfo(t_as_mat)
end

function split_long_edges(x::Vector, mi::MeshInfo, l::Float64)
	#= Split edges longer than a factor greater than the mean length

	INPUT
	x - vector of positions
	mi - the mesh
	l - splits edges longer than l

	OUTPUT
	new_x - new vector of positions
	new_mi - new MeshInfo object
	=#

	new_x = copy(x)
	new_t = reshape(copy(mi.t), length(mi.t))

	# maps edges to triangles to which they belong
	edge_map = construct_edge_map(mi.t)

	# queue of triangles to visit
	q = [i for i in 1:size(mi.t,2)]

	function get_other_point_in_tri(tind, edge)
		if new_t[3*tind-2] != edge[1] && new_t[3*tind-2] != edge[2]
			return new_t[3*tind-2]
		elseif new_t[3*tind-1] != edge[1] && new_t[3*tind-1] != edge[2]
			return new_t[3*tind-1]
		else
			return new_t[3*tind]
		end
	end

	function replace_tri(edge, oldtri, newtri)
		(t1,t2) = edge_map[edge]

		push!(edge_map, edge => (t1 == oldtri ? (newtri, t2) : (t1, newtri)))
	end

	i2e = (a,b) -> (min(a,b), max(a,b))

	while !isempty(q)
		tind = popfirst!(q)

		i1 = new_t[3*tind-2]; i2 = new_t[3*tind-1]; i3 = new_t[3*tind];

		# find longest edge to split
		l1 = norm(new_x[3*i2-2:3*i2] .- new_x[3*i1-2:3*i1])
		l2 = norm(new_x[3*i3-2:3*i3] .- new_x[3*i2-2:3*i2])
		l3 = norm(new_x[3*i1-2:3*i1] .- new_x[3*i3-2:3*i3])

		edge_length = max(l1, l2, l3)
		e_split = [(i1, i2), (i2, i3), (i3, i1)][argmax([l1, l2, l3])]
		e_split = i2e(e_split...)
		(i1, i2) = e_split												# i1, i2 overwritten here

		# split edge if too long
		if edge_length > l
			(t1, t2) = edge_map[e_split]

			# remove edge from edge map
			delete!(edge_map, e_split)

			# create new point
			append!(new_x, (new_x[3*i2-2:3*i2] .+ new_x[3*i1-2:3*i1])/2)
			new_i = length(new_x) ÷ 3

			# overwrite triangles and update edge_map
			t1_p = get_other_point_in_tri(t1, e_split)			# point in t1 not included in split edge
			new_t[3*t1-2:3*t1] .= [i1, t1_p, new_i]	
			append!(new_t, [i2, t1_p, new_i])
			t1_add = length(new_t) ÷ 3					# index of newly created triangle	
			replace_tri(i2e(i1, t1_p), t1, t1)
			replace_tri(i2e(i2, t1_p), t1, t1_add)
			push!(edge_map, i2e(t1_p, new_i) => (t1, t1_add))

			if t2 > 0
				t2_p = get_other_point_in_tri(t2, e_split)			# point in t1 not included in split edge
				new_t[3*t2-2:3*t2] .= [i1, t2_p, new_i]	
				append!(new_t, [i2, t2_p, new_i])
				t2_add = length(new_t) ÷ 3					# index of newly created triangle	
				replace_tri(i2e(i1, t2_p), t2, t2)
				replace_tri(i2e(i2, t2_p), t2, t2_add)
				push!(edge_map, i2e(t2_p, new_i) => (t2, t2_add))

				# add middle edges
				push!(edge_map, i2e(i1, new_i) => (t1, t2))
				push!(edge_map, i2e(i2, new_i) => (t1_add, t2_add))
			else
				# add middle edges (with no second triangle)
				push!(edge_map, i2e(i1, new_i) => (t1, 0))
				push!(edge_map, i2e(i2, new_i) => (t1_add, 0))
			end

			# add new triangles to queue (this leads to some duplicates, but that's fine)
			push!(q, t1, t1_add)
			t2 > 0 && push!(q, t2, t2_add)
		end
	end

	# construct new mesh
	t_as_mat = reshape(new_t, (3, length(new_t) ÷ 3))

	tnbor = triangle_connectivity(t_as_mat)
	order_triangles!(t_as_mat, tnbor)

	return new_x, MeshInfo(t_as_mat)
end

function edge_flips(x::Vector{Float64}, mi::MeshInfo)
	#= flip edges to improve valence of each point (optimal is 6 for interior and 4 for boundary)
	
	INPUT
	mi - MeshInfo object

	OUTPUT
	new_mi - a new MeshInfo object
	=#

	nbor_count = copy(mi.nbor_count)
	t = copy(mi.t)
	isbndry = mi.isbndry

	tnbor = copy(mi.tnbor)

	# maps edges to triangles to which they belong
	edge_map = construct_edge_map(mi.t)

	# return triangles neighboring t1 that aren't t2
	@inline function get_other_tris(t1, t2)
		t11, t12, t13 = tnbor[1:3, t1]
		t11 == t2 && return t12, t13
		t12 == t2 && return t11, t13
		return t11, t12
	end

	# replace neighboring triangle t2 of t1 with t3
	@inline function replace_tri_in_tnbor(t1, t2, t3)
		if t1 == 0
			return
		end

		@assert t2 != 0
		@assert t2 ∈ tnbor[:,t1]

		tnbor[1,t1] == t2 && (tnbor[1,t1] = t3)
		tnbor[2,t1] == t2 && (tnbor[2,t1] = t3)
		tnbor[3,t1] == t2 && (tnbor[3,t1] = t3)
	end

	function get_other_point_in_tri(tind, edge)
		if t[1,tind] != edge[1] && t[1,tind] != edge[2]
			return t[1,tind]
		elseif t[2,tind] != edge[1] && t[2,tind] != edge[2]
			return t[2,tind]
		else
			return t[3,tind]
		end
	end

	function replace_tri(edge, oldtri, newtri)
		(t1,t2) = edge_map[edge]

		push!(edge_map, edge => (t1 == oldtri ? (newtri, t2) : (t1, newtri)))
	end

	i2e = (a,b) -> (min(a,b), max(a,b))

	# iterate through all edges
	for (edge, (t1, t2)) in edge_map
		(i1, i2) = edge

		# only look at edges that between two triangles
		if t2 == 0
			continue
		end

		(i1, i2) = edge
		t1_p = get_other_point_in_tri(t1, edge)
		t2_p = get_other_point_in_tri(t2, edge)

		cur_valence = (nbor_count[i1] - 6 + 2*isbndry[i1])^2 + (nbor_count[i2] - 6 + 2*isbndry[i2])^2
		cur_valence += (nbor_count[t1_p] - 6 + 2*isbndry[t1_p])^2 + (nbor_count[t2_p] - 6 + 2*isbndry[t2_p])^2

		new_valence = (nbor_count[i1]-1 - 6 + 2*isbndry[i1])^2 + (nbor_count[i2]-1 - 6 + 2*isbndry[i2])^2
		new_valence += (nbor_count[t1_p]+1 - 6 + 2*isbndry[t1_p])^2 + (nbor_count[t2_p]+1 - 6 + 2*isbndry[t2_p])^2

		t11, t12 = get_other_tris(t1, t2)
		t21, t22 = get_other_tris(t2, t1)

		# ensure that the neighboring triangles don't equal each other
		if t11 == t21 || t12 == t22 || t11 == t12 || t21 == t22 || t21 == t12 || t11 == t22
			continue
		end

		# flip edges if new valence is better
		if new_valence < cur_valence
			nbor_count[i1] -= 1
			nbor_count[i2] -= 1
			nbor_count[t1_p] += 1
			nbor_count[t2_p] += 1

			# update triangles
			t[:,t1] .= [i1, t1_p, t2_p]
			t[:,t2] .= [i2, t1_p, t2_p]

			# ------------ udpate tnbor ----------------

			# swap definitions so that t11 and t21 are the triangles touching i1
			if i1 ∉ t[:,t11]
				swp = t11;
				t11 = t12;
				t12 = swp;
			end

			if i1 ∉ t[:,t21]
				swp = t21;
				t21 = t22;
				t22 = swp;
			end

			replace_tri_in_tnbor(t21, t2, t1)
			replace_tri_in_tnbor(t12, t1, t2)
			tnbor[:,t1] .= [t2, max(t11, t21), min(t11, t21)]	# ensures 0 is last if present
			tnbor[:,t2] .= [t1, max(t12, t22), min(t12, t22)]


			# update edge map with the new edge and the other two that have changed
			delete!(edge_map, edge)
			push!(edge_map, i2e(t1_p, t2_p) => (t1, t2))
			replace_tri(i2e(i1, t2_p), t2, t1)
			replace_tri(i2e(i2, t1_p), t1, t2)
		end
	end

	try
		tnbor = triangle_connectivity(t)
		order_triangles!(t, tnbor)
	catch err
		@assert false "Mesh degenerated"
	end

	return MeshInfo(t)
end

function vertex_shift(x::Vector{Float64}, mi::MeshInfo, r::Float64=0.4)
	#= shift the vertices to the center of its 1-neighbor ring centroid

	INPUT
	x - vertices
	mi - MeshInfo object
	r - fraction of the way to adjust vertices to centroid of neighboring vertices

	OUTPUT
	x - new vertices
	=#

	xc = copy(x)

	for i1 = 1:length(x) ÷ 3
		cent = zeros(3)

		# skip if boundary point
		if mi.isbndry[i1]
			continue
		end

		for j = 1:mi.nbor_count[i1]
			i2 = mi.nbor[j,i1]
			cent .+= x[3*i2-2:3*i2]
		end

		cent ./= mi.nbor_count[i1]

		xc[3*i1-2:3*i1] .= (1-r).*xc[3*i1-2:3*i1] .+ r.*cent
	end

	return xc
end

function remesh(x::Vector{Float64}, mi::MeshInfo, l_lo, l_hi::Float64, projection_x::Vector{Float64} = x, projection_mi::MeshInfo = mi)
	#= perform remeshing by iteratively collapsing short edges and splitting long edges

	INPUT
	x - positions of mesh points
	mi - MeshInfo object
	l_lo - length below which edges are collapsed
	l_hi - length above which edges are split
	projection_x - optional additional set of vertices used to perform projection procedure
	projection_mi - MeshInfo object to go along with projection_x

	OUTPUT
	x_new - new set of vertices
	mi_new - new MeshInfo object
	=#

	prev_size_t = (3, 0)
	prev_length_x = 0

	# build octree for projection
	o = create_octree_from_mesh(projection_x, projection_mi);

	x_new = copy(x)
	mi_new = MeshInfo(mi.t)

	# iterate until no more changes
	for i = 1:10
		x_new, mi_new = split_long_edges(x_new, mi_new, l_hi)

		x_new, mi_new = collapse_short_edges(x_new, mi_new, l_lo)

		mi_new = edge_flips(x_new, mi_new)

		x_new = vertex_shift(x_new, mi_new, i <= 5 ? 0.4 : 1.0) # adjust shifting parameter

		# project points onto original mesh
		for j = 1:length(x_new)÷3
			new_p, _ = find_closest_point(x_new[3*j-2:3*j], o, projection_x, projection_mi)
			x_new[3*j-2:3*j] .= new_p
		end
	end

	return x_new, mi_new
end

# ====================================== FILE WRITING AND READING ===========================================




function save_mesh(mi::MeshInfo, filename::String)
	#= Save the mesh data to a file. Only writes the triangle information since
	MeshInfo objects can be completely reconstructed from this alone.

	INPUT
	mi - MeshInfo object
	filename - the filename

	=#

	open(filename, "w") do io
		# write version number
		version = 1
		write(io, Int64(version))

		# write sizes of matrices
		Nt = size(mi.t, 2)
		write(io, Int64(Nt))

		# write data to file
		for i in eachindex(mi.t)
			write(io, mi.t[i])
		end
	end
end

function read_mesh(filename::String)
	#= Reads in a MeshInfo object from a file.

	INPUT
	filename - the filename

	OUTPUT
	a new MeshInfo object

	=#

	open(filename, "r") do io
		# write version number
		version = read(io, Int64)

		Nt = read(io, Int64)

		# read in triangle information
		t = Array{Int64}(undef, 3, Nt)
		read!(io, t)

		return MeshInfo(t)
	end
end

function save_vector(x::Vector{T}, filename::String) where T
	#= Saves a vector to a file

	INPUT
	x - MeshInfo object
	filename - the filename

	=#

	open(filename, "w") do io
		# write size of vector
		write(io, Int64(length(x)))

		# write data to file
		for i in eachindex(x)
			write(io, x[i])
		end
	end
end

function read_vector(filename::String, ::Type{T}) where T
	#= Reads in a vector of type T from a file

	INPUT
	filename - the filename

	=#

	open(filename, "r") do io
		# write size of vector
		sz = read(io, Int64)

		x = Array{T}(undef, sz)

		# read in data
		read!(io, x)

		return x
	end
end
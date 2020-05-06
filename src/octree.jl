struct Rect{T}
	x::Tuple{Float64,Float64,Float64}
	w::Tuple{Float64,Float64,Float64}
	el::T
end

Rect(x,w) = Rect(x,w,nothing)

function intersects(r1::Rect, r2::Rect)
	return all(r1.x .<= (r2.x .+ r2.w)) && all(r2.x .<= (r1.x .+ r1.w))
end

# ------------------- OctNode ------------------

mutable struct OctNode{T}
	xyz::Union{OctNode,Nothing}
	xyz_::Union{OctNode,Nothing}
	xy_z::Union{OctNode,Nothing}
	xy_z_::Union{OctNode,Nothing}
	x_yz::Union{OctNode,Nothing}
	x_yz_::Union{OctNode,Nothing}
	x_y_z::Union{OctNode,Nothing}
	x_y_z_::Union{OctNode,Nothing}
	c::Tuple{Float64,Float64,Float64}	# center
	w::Tuple{Float64,Float64,Float64}	# width
	el::Vector{Rect{T}}
end

MeshOctNode = OctNode{Int64}

function OctNode{T}() where T
	return OctNode([nothing for i=1:8]..., (0.0,0.0,0.0), (0.0,0.0,0.0), Rect{T}[])
end

function OctNode{T}(c::Tuple{Float64,Float64,Float64}, w::Tuple{Float64,Float64,Float64}) where T
	return OctNode([nothing for i=1:8]..., c, w, Rect{T}[])
end

function children(o::OctNode{T})::Array{OctNode{T},1} where T
	if o.xyz == nothing
		return OctNode{T}[]
	else
		return [o.xyz, o.xyz_, o.xy_z, o.xy_z_, o.x_yz, o.x_yz_, o.x_y_z, o.x_y_z_]
	end
end

has_children(o::OctNode) = o.xyz != nothing;

depth(o::OctNode, n=0) = has_children(o) ? maximum([depth(i, n+1) for i in children(o)]) : n;

# constructor for octnode
Rect(o::OctNode) = Rect(o.c .- o.w./2, o.w);

Base.in(r::Rect, o::OctNode) = intersects(r, Rect(o));

# iterate through all rects of children
function Base.iterate(o::OctNode{T}, state=([o], 1)) where T
	q, i = state

	if isempty(q)
		return nothing
	end

	node = first(q)

	while !isempty(q)
		node = first(q)

		if has_children(node)
			popfirst!(q)
			append!(q, children(node))
		else
			break
		end
	end

	if i <= length(node.el)
		return (node.el[i], (q, i+1))
	else
		popfirst!(q)
		return iterate(o, (q, 1))
	end
end

Base.length(o::OctNode) = has_children(o) ? sum([length(i) for i in children(o)]) : length(o.el);
Base.eltype(o::OctNode{T}) where T = Rect{T} 

# --------------- main functions to insert and split nodes -------------------------

function split!(o::OctNode{T}, split_max::Int64) where T
	if length(o.el) > split_max
		o.xyz    = OctNode{T}(o.c .+ o.w.*( 1, 1, 1)./4, o.w./2)
		o.xyz_   = OctNode{T}(o.c .+ o.w.*( 1, 1,-1)./4, o.w./2)
		o.xy_z   = OctNode{T}(o.c .+ o.w.*( 1,-1, 1)./4, o.w./2)
		o.xy_z_  = OctNode{T}(o.c .+ o.w.*( 1,-1,-1)./4, o.w./2)
		o.x_yz   = OctNode{T}(o.c .+ o.w.*(-1, 1, 1)./4, o.w./2)
		o.x_yz_  = OctNode{T}(o.c .+ o.w.*(-1, 1,-1)./4, o.w./2)
		o.x_y_z  = OctNode{T}(o.c .+ o.w.*(-1,-1, 1)./4, o.w./2)
		o.x_y_z_ = OctNode{T}(o.c .+ o.w.*(-1,-1,-1)./4, o.w./2)

		# insert rectangles into children nodes now
		while !isempty(o.el)
			insert!(o, pop!(o.el), split_max, false)
		end
	end
end

function insert!(o::OctNode, r::Rect, split_max::Int64, should_split=true)
	if has_children(o)
		for child in children(o)
			insert!(child, r, split_max)
		end
	else
		if r ∈ o
			push!(o.el, r)
			should_split && split!(o, split_max)
		end
	end
end

# ------------------- triangle utilities ------------------------------------------

function closest_point_on_triangle(p_3d, v1_3d, v2_3d, v3_3d)
	# returns closest point and distance

	v1 = (v1_3d .- v3_3d);
	v2 = (v2_3d .- v3_3d);
	p_tmp = (p_3d .- v3_3d);

	# calculate 3d normal
	n_3d = cross(v1, v2)
	n_3d ./= norm(n_3d)

	# project p
	p = p_tmp .- dot(p_tmp, n_3d).*n_3d

	# matrix for projections
	A = [v1 v2]

	# in triangle
	b = A \ p
	if all(b .> 0.0) && sum(b) <= 1
		proj = A*b .+ v3_3d
		return proj, norm(proj .- p_3d)
	end


	# edge 1
	b = v1 \ p;
	proj1 = v1*b .+ v3_3d
	d1 = Inf

	if b > 0.0 && b <= 1
		d1 = norm(proj1 .- p_3d)
	end

	# edge 2
	b = v2 \ p;
	proj2 = v2*b .+ v3_3d
	d2 = Inf

	if b > 0.0 && b <= 1
		d2 = norm(proj2 .- p_3d)
	end

	# edge 3
	b = (v2 .- v1) \ (p .- v1)
	proj3 = (v2 .- v1)*b .+ v1 .+ v3_3d
	d3 = Inf

	if b > 0.0 && b <= 1
		d3 = norm(proj3 .- p_3d)
	end

	# one of the vertices
	dp1 = norm(p_3d .- v1_3d)
	dp2 = norm(p_3d .- v2_3d)
	dp3 = norm(p_3d .- v3_3d)

	ind = argmin([d1, d2, d3, dp1, dp2, dp3])
	if ind == 1
		return proj1, d1
	elseif ind == 2
		return proj2, d2
	elseif ind == 3
		return proj3, d3
	elseif ind == 4
		return v1_3d, dp1
	elseif ind == 5
		return v2_3d, dp2
	else
		return v3_3d, dp3
	end
end

function find_closest_point(p::Vector{Float64}, o::MeshOctNode, x::Vector{Float64}, mi::MeshInfo)
	#= find the closest point on a mesh by searching through the octree

	INPUT
	p - input point
	o - octree
	x - positions of nodes
	mi - MeshInfo object

	OUTPUT
	closest - the closest point to p that lies on the mesh
	d - distance to the point
	=#

	# find cell in which p resides
	p_rect = Rect(Tuple(p), (0.,0.,0.))

	if p_rect ∉ o
		throw(DomainError(p, "Point outside of octree"))
	end

	close_node = o
	while has_children(close_node)
		found_rects = false;
		for child in children(close_node)
			if p_rect ∈ child && length(child) > 0
				close_node = child
				found_rects = true
				break
			end
		end

		if !found_rects
			break
		end
	end

	# find closest triangle
	d = Inf
	closest = zeros(3)

	for r in close_node
		tind = r.el
		i1 = mi.t[1,tind]; i2 = mi.t[2,tind]; i3 = mi.t[3,tind]; 
		closest_tmp, d_tmp = closest_point_on_triangle(p, x[3*i1-2:3*i1], x[3*i2-2:3*i2], x[3*i3-2:3*i3])

		if d_tmp < d
			d = d_tmp
			closest .= closest_tmp
		end
	end

	# now search all triangles in cells that are within d of point p
	p_rect = Rect(Tuple(p) .- d, (2*d,2*d,2*d))

	q = [o]

	while !isempty(q)
		node = popfirst!(q)

		if has_children(node)
			# append children to queue if p_rect touches them
			for child in children(node)
				if p_rect ∈ child
					push!(q, child)
				end
			end
		elseif node !== close_node
			# search through triangles of these nodes
			for r in node
				if !intersects(p_rect, r)
					# ignore triangles whose bounding boxes don't overlap p_rect
					continue
				end

				tind = r.el
				i1 = mi.t[1,tind]; i2 = mi.t[2,tind]; i3 = mi.t[3,tind]; 
				closest_tmp, d_tmp = closest_point_on_triangle(p, x[3*i1-2:3*i1], x[3*i2-2:3*i2], x[3*i3-2:3*i3])

				if d_tmp < d
					d = d_tmp
					closest .= closest_tmp
				end
			end

		end
	end

	return closest, d
end

function create_octree_from_mesh(x::Vector{Float64}, mi::MeshInfo)
	#= create an octree
	divide until each child node intersects with less than 2 x max(neighbor_count) triangle
	bounding boxes

	INPUT
	x - positions of nodes
	mi - MeshInfo

	OUTPUT
	octtree - an octree of type MeshOctNode
	=#

	maxN = 2*maximum(mi.nbor_count)

	# calculate bounding box
	box_min_x = reduce(min, x[i] for i in 1:3:length(x));
	box_min_y = reduce(min, x[i] for i in 2:3:length(x));
	box_min_z = reduce(min, x[i] for i in 3:3:length(x));
	box_max_x = reduce(max, x[i] for i in 1:3:length(x));
	box_max_y = reduce(max, x[i] for i in 2:3:length(x));
	box_max_z = reduce(max, x[i] for i in 3:3:length(x));

	box_c = ((box_min_x + box_max_x)/2, (box_min_y + box_max_y)/2, (box_min_z + box_max_z)/2)
	box_w = (1.5*(box_max_x - box_min_x), 1.5*(box_max_y - box_min_y), 1.5*(box_max_z - box_min_z)) # add 50% to be safe

	function rect_from_tri(ind)
		i1 = mi.t[1,ind]; i2 = mi.t[2,ind]; i3 = mi.t[3,ind];
		x1 = x[3*i1-2]; y1 = x[3*i1-1]; z1 = x[3*i1];
		x2 = x[3*i2-2]; y2 = x[3*i2-1]; z2 = x[3*i2];
		x3 = x[3*i3-2]; y3 = x[3*i3-1]; z3 = x[3*i3];

		min_x = min(x1, x2, x3); min_y = min(y1, y2, y3); min_z = min(z1, z2, z3);
		max_x = max(x1, x2, x3); max_y = max(y1, y2, y3); max_z = max(z1, z2, z3);

		return Rect((min_x, min_y, min_z), (max_x, max_y, max_z) .- (min_x, min_y, min_z), ind)
	end

	root = MeshOctNode(box_c, box_w)

	for i = 1:size(mi.t, 2)
		insert!(root, rect_from_tri(i), maxN)
	end

	return root
end
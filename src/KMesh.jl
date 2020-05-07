module KMesh

using Requires

struct MeshInfo
	t::Array{Int64,2}
	nbor::Array{Int64,2}
	nbor_count::Array{Int64,1}
	tmem::Array{Int64,2}
	tmem_count::Array{Int64,1}
	isbndry::BitArray{1}
	bndry_path::Array{Int64,1}
	tnbor::Array{Int64,2}
end

function MeshInfo(t)
	Np = maximum(t)

	nbor, nbor_count, tmem, tmem_count = calculate_neighbors(t, Np)
	isbndry = find_boundary_vertices(nbor_count, tmem_count);
	bndry_path = find_boundary_edges(isbndry, t);
	tnbor = triangle_connectivity(t)

	return MeshInfo(t, nbor, nbor_count, tmem, tmem_count, isbndry, bndry_path, tnbor)
end

include("octree.jl")
include("mesh.jl")
include("curvature.jl")

# export types
export MeshInfo, Rect, OctNode

# mesh.jl functions
export order_triangles!, calculate_neighbors, triangle_connectivity, find_boundary_vertices,
	find_boundary_edges, contour_length, area_volume, area, areas, normals, subdivide,
	mean_edge_length, construct_edge_map, collapse_short_edges, split_long_edges, edge_flips, vertex_shift,
	remesh, save_mesh, read_mesh, save_vector, read_vector

# curvature.jl functions
export calculate_curvature, calculate_curvature!

# octree.jl functions
export find_closest_point, create_octree_from_mesh

# include plotting code
function __init__()
	@require MATLAB="10e44e05-a98a-55b3-a45b-ba969058deb6" include("plotting.jl")
end

end # module

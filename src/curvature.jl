using LinearAlgebra
# import ReverseDiff

# import LinearAlgebra.normalize
# import LinearAlgebra.normalize!

# # define normalization functions for vectors of duals
# @inline function normalize(v)
# 	return copy(v) ./ norm(v);
# end

# @inline function normalize!(v)
# 	v ./= norm(v);
# end

function calculate_curvature!(p, t, nbor, nbor_count, tmem, tmem_count, isbndry, H, K, N, A)
	#= Calculates the mean curvature at every point as well as normal vectors
	Overwrites the storage H, K, N, and A
	
	INPUT
	p - a stacked vector of coordinates [x1 y1 z1 x2 y2 z2 ...]
	t - indices of triangles
	nbor - list of neighbors
	nbor_count - vector of neighbor numbers
	tmem - triangle member matrix
	tmem_count - number of triangles to which each vertex belongs
	isbndry - vector of booleans if a vertex is on the border
	
	OUTPUT
	H - mean curvatures at each vertex
	K - Gaussian curvature at each vertex
	N - 3 x Np matrix of normal vectors at each vertex
	A - areas of each vertex
	=#
	
	Np = length(p) ÷ 3

	# fill with zeros
	fill!(N, zero(eltype(p)))
	fill!(A, zero(eltype(p)))
	fill!(K, one(eltype(p))*π)
	K .*= (2 .- isbndry)
	
	Atri_tot = 0.0

	v12 = similar(p, 3)
	v13 = similar(p, 3)
	v23 = similar(p, 3)
	
	# iterate through each triangle
	for i = 1:size(t, 2)
		
		i1 = t[1,i]; i2 = t[2,i]; i3 = t[3,i];
		
		# calculate bond vectors
		v12 .= p[3*i2-2:3*i2] .- p[3*i1-2:3*i1]
		v13 .= p[3*i3-2:3*i3] .- p[3*i1-2:3*i1]
		v23 .= p[3*i3-2:3*i3] .- p[3*i2-2:3*i2]
		
		nv12_2 = dot(v12, v12)
		nv13_2 = dot(v13, v13)
		nv23_2 = dot(v23, v23)
		
		# dot products of bond vectors
		d1 =  dot(v12, v13) / sqrt(nv12_2*nv13_2)
		d2 = -dot(v12, v23) / sqrt(nv12_2*nv23_2)# negative due to direction
		d3 =  dot(v13, v23) / sqrt(nv13_2*nv23_2)

		K[i1] -= acos(d1)
		K[i2] -= acos(d2)
		K[i3] -= acos(d3)
		
		# cotangents of angles
		cot1 = d1 / sqrt(1 - d1^2)
		cot2 = d2 / sqrt(1 - d2^2)
		cot3 = d3 / sqrt(1 - d3^2)
		
		# determine if triangle is obtuse
		# obtuse = (d1 <= 0.0 || d2 <= 0.0 || d3 <= 0.0)
		# obtuse = 1.0 - sign(max(d1*d2*d3, 0.0))

		Atri = norm(cross(v12, v13))/2 # area of the triangle
		Atri_tot += Atri
		
		# if obtuse
		# 	# add 1/2 area if angle is obtuse, otherwise add 1/4 area
		# 	A[i1] += (d1 <= 0.0 ? Atri/2 : Atri/4)
		# 	A[i2] += (d2 <= 0.0 ? Atri/2 : Atri/4)
		# 	A[i3] += (d3 <= 0.0 ? Atri/2 : Atri/4)
		# else
		# 	# triangle not obtuse so use Voronoi areas
		# 	A[i1] += cot2/8 * nv13_2
		# 	A[i1] += cot3/8 * nv12_2
			
		# 	A[i2] += cot1/8 * nv23_2
		# 	A[i2] += cot3/8 * nv12_2
			
		# 	A[i3] += cot1/8 * nv23_2
		# 	A[i3] += cot2/8 * nv13_2
		# end

		# a1, a2, a3 = assign_area(d1, d2, d3, Atri, cot1, cot2, cot3, nv12_2, nv13_2, nv23_2)
		# A[i1] += a1
		# A[i2] += a2
		# A[i3] += a3

		# a1 = (cot2*nv13_2 + cot3*nv12_2)/8
		# a2 = (cot1*nv23_2 + cot3*nv12_2)/8
		# a3 = (cot1*nv23_2 + cot2*nv13_2)/8

		# ob_a1 = Atri/3
		# ob_a2 = Atri/3
		# ob_a3 = Atri/3

		# obtuse = (1.0 - tanh(20*d1*d2*d3))/2
		# A[i1] += obtuse*ob_a1 + (1.0 - obtuse)*a1
		# A[i2] += obtuse*ob_a2 + (1.0 - obtuse)*a2
		# A[i3] += obtuse*ob_a3 + (1.0 - obtuse)*a3
		A[i1] += Atri/3
		A[i2] += Atri/3
		A[i3] += Atri/3

		# obtuse = aux1(d1) + aux1(d2) + aux1(d3)
		# A[i1] += obtuse*aux2(d1, Atri) + (1.0 - obtuse)*a1
		# A[i2] += obtuse*aux2(d2, Atri) + (1.0 - obtuse)*a2
		# A[i3] += obtuse*aux2(d3, Atri) + (1.0 - obtuse)*a3


		# add contribution to normal vectors (negatives ensure direction is correct)
		N[:,i1] .+= -cot2 * v13
		N[:,i1] .+= -cot3 * v12

		N[:,i2] .+= -cot1 * v23
		N[:,i2] .+= cot3 * v12

		N[:,i3] .+= cot1 * v23
		N[:,i3] .+= cot2 * v13
		
	end
	
	# divide Gaussian curvature by the mixed areas
	K ./= A
	
	new_n = zeros(eltype(p), 3)

	# calculate mean curvature and normalize normals
	# find any normals with 0 mean curvature and simply take the normal of a triangle touching it
	for i = 1:Np
		this_norm = sqrt(N[1,i]^2 + N[2,i]^2 + N[3,i]^2 + eps()^2)

		if this_norm <= sqrt(2)*eps()
			t_ind = tmem[1,i] # take the first triangle arbitrarily
			
			i1 = t[1,t_ind]; i2 = t[2,t_ind]; i3 = t[3,t_ind]
			# calculate bond vectors
			v12 .= p[3*i2-2:3*i2] .- p[3*i1-2:3*i1]
			v13 .= p[3*i3-2:3*i3] .- p[3*i1-2:3*i1]

			new_n .= cross(v12, v13);
		
			N[:,i] .= new_n
			N[:,i] ./= norm(new_n)
		else
			N[:,i] ./= this_norm
		end

		H[i] = this_norm / (4*A[i])
	end
	
	# fix curvature and normals for the border vertices
	# use Max method to estimate vertex normals
	p_i = zeros(eltype(p), 3)

	for i = 1:Np
		if !isbndry[i]
			continue
		end
		
		fill!(new_n, zero(eltype(p)))
		p_i = p[3*i-2:3*i]
		
		for j = 1:tmem_count[i]
			tind = tmem[j,i]
			
			i1 = t[1,tind]; i2 = t[2,tind]; i3 = t[3,tind];

			# rearrange vertices so i1 is i while maintaining consistent normals
			if i2 == i
				i2 = i3; i3 = i1;
			elseif i3 == i
				i3 = i2; i2 = i1;
			end
		
			# calculate bond vectors
			v12 .= p[3*i2-2:3*i2] .- p_i
			v13 .= p[3*i3-2:3*i3] .- p_i
			
			new_n .+= cross(v12, v13) / (dot(v12,v12)*dot(v13,v13))
		end
		
		new_n ./= norm(new_n)
		
		# change direction to match that of an interior neighbor
		j = 1
		while isbndry[nbor[j,i]]
			j += 1
		end
		
		if dot(new_n, N[:,nbor[j,i]]) < 0
			new_n .*= -1
		end
		
		N[:,i] .= new_n
		
		# # calculate curvature (using finite differences projection)
		# # h_new = 1.0
		# h_new = 0.0
		
		# for j = 1:nbor_count[i]
		# 	nind = nbor[j,i]
		# 	diff_p = p[3*nind-2:3*nind] .- p_i
		# 	norm_diff_p_2 = dot(diff_p, diff_p)
			
		# 	sgn = 1;
		# 	if dot(N[:,nind], new_n) < 0
		# 		sgn = -1
		# 	end
			
		# 	# h_new *= abs(dot(sgn*N[:,nind] .- new_n, diff_p)) / norm_diff_p_2
		# 	h_new += abs(dot(sgn*N[:,nind] .- new_n, diff_p)) / norm_diff_p_2
		# end

		h_new = zero(eltype(p))
		weights = zero(eltype(p))

		for j = 1:nbor_count[i]
			nind = nbor[j,i]

			if isbndry[nind]
				continue
			end

			diff_p = p[3*nind-2:3*nind] .- p_i
			norm_diff_p_2 = dot(diff_p, diff_p)
			
			h_new += H[nind]/norm_diff_p_2
			weights += 1/norm_diff_p_2
		end
		
		# find arithmetic mean of the curvatures
		# h_new = h_new^(1/nbor_count[i])
		H[i] = h_new / weights
	end

	return 0
end

# ReverseDiff.@forward function aux1(a)
# 	# q = 0.1

# 	# if a < 0
# 	# 	return 1.0
# 	# elseif a < q
# 	# 	return ((a - q)/q)^2
# 	# end

# 	# return 0.0
# 	return (a <= 0 ? 1.0 : 0.0)
# 	# return (1 - tanh(20*a))/2
# end

# ReverseDiff.@forward function aux2(di, Atri)
# 	return (di <= 0.0 ? Atri/2 : Atri/4)
# 	# return Atri/3
# end

@inline function calculate_curvature(p, mi::MeshInfo)
	H = zeros(eltype(p), length(p)÷3)
	K = zeros(eltype(p), length(p)÷3)
	N = zeros(eltype(p), 3, length(p)÷3)
	A = zeros(eltype(p), length(p)÷3)

	calculate_curvature!(p, mi.t, mi.nbor, mi.nbor_count, mi.tmem, mi.tmem_count, mi.isbndry, H, K, N, A)

	return H, K, N, A
end

@inline function calculate_curvature!(p, mi::MeshInfo, H, K, N, A)
	calculate_curvature!(p, mi.t, mi.nbor, mi.nbor_count, mi.tmem, mi.tmem_count, mi.isbndry,
	H, K, N, A)
end
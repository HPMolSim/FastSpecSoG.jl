# precomputings

@inline @inbounds function set_zeros!(a::AbstractArray{T}) where{T <: Number}
    return fill!(a, zero(T))
end

@inline function chebpoly(n::Int, x::T, scale::T) where{T}
    return cos(n * acos(x / scale))
end

# precompute the parameters to be used
function FFCT_precompute(L::NTuple{3, T}, N_grid::NTuple{3, Int}, uspara::USeriesPara{T}, M_mid::Int, n_atoms::Int) where{T}
    L_x, L_y, L_z = L
    N_x, N_y, N_z = N_grid
    k_x = Vector{T}([2π * i / L_x for i in -N_x:N_x])
    k_y = Vector{T}([2π * i / L_y for i in -N_y:N_y])
    r_z = Vector{T}([L_z / 2 * ( 1.0 - cos((2i - 1)*π / 2N_z) ) for i in 1:N_z])
    
    k_mat = zeros(T, 2N_x + 1, 2N_y + 1)
    kx_mat = zeros(T, 2N_x + 1, 2N_y + 1)
    ky_mat = zeros(T, 2N_x + 1, 2N_y + 1)

    for i in 1:2N_x + 1, j in 1:2N_y + 1
        k_mat[i, j] = k_x[i]^2 + k_y[j]^2
        kx_mat[i, j] = k_x[i]
        ky_mat[i, j] = k_y[j]
    end

    # xy k space grid, z real space grid
    H_r = zeros(Complex{T}, 2N_x + 1, 2N_y + 1, N_z)
    # xy k space grid, z ChebCoefs
    H_c = zeros(Complex{T}, 2N_x + 1, 2N_y + 1, N_z)
    # solution to the linear eqs, xy k space grid, z ChebCoefs
    H_s = zeros(Complex{T}, 2N_x + 1, 2N_y + 1, N_z)

    M = length(uspara.sw)
    M_range = M - M_mid
    us_mat = zeros(T, 2N_x + 1, 2N_y + 1, M_range)
    for l in 1:M_range
        sl, wl = uspara.sw[l + M_mid]
        us_mat[:, :, l] = exp.( - sl^2 .* k_mat ./ 4)
    end

    ivsm = inverse_mat(N_grid, L[3], k_x, k_y)

    # the boundary condition at z = 0 and z = L_z
    b_l = zeros(Complex{T}, 2N_x + 1, 2N_y + 1)
    b_u = zeros(Complex{T}, 2N_x + 1, 2N_y + 1)
    phase_x = zeros(Complex{T}, 2N_x + 1)
    phase_y = zeros(Complex{T}, 2N_y + 1)

    phase_xs = zeros(Complex{T}, 2N_x + 1, n_atoms)
    phase_ys = zeros(Complex{T}, 2N_y + 1, n_atoms)
    phase_xys = zeros(Complex{T}, 2N_x + 1, 2N_y + 1, n_atoms)

    rhs = zeros(Complex{T}, N_z + 2)
    sol = zeros(Complex{T}, N_z + 2)

    sort_z = zeros(Int, n_atoms)
    z = zeros(T, n_atoms)

    size_dict = Dict('i' => 2N_x + 1, 'j' => 2N_y + 1, 'k' => N_z, 'l' => M_range, 'n' => n_atoms)
    temp_ijlk = zeros(Complex{T}, (size_dict['i'], size_dict['j'], size_dict['l'], size_dict['k']))
    temp_ijl = zeros(Complex{T}, (size_dict['i'], size_dict['j'], size_dict['l']))
    z_coef = zeros(T, M_range, N_z, n_atoms)
    exp_coef = zeros(T, M_range, N_z, n_atoms)
    
    return k_x, k_y, k_mat, r_z, Complex{T}.(us_mat), H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef
end

@inbounds function revise_phase_neg!(phase_x::AbstractArray{Complex{T}, 1}, phase_y::AbstractArray{Complex{T}, 1}, k_x::Vector{T}, k_y::Vector{T}, x::T, y::T) where{T}

    phase_x .= exp.(-T(1)im .* k_x .* x)
    phase_y .= exp.(-T(1)im .* k_y .* y)

    return nothing
end

@inbounds function revise_phase_pos!(phase_x::AbstractArray{Complex{T}, 1}, phase_y::AbstractArray{Complex{T}, 1}, k_x::Vector{T}, k_y::Vector{T}, x::T, y::T) where{T}

    phase_x .= exp.(T(1)im .* k_x .* x)
    phase_y .= exp.(T(1)im .* k_y .* y)

    return nothing
end

@inbounds function revise_phase_neg_all!(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, phase_xs::Array{Complex{T}, 2}, phase_ys::Array{Complex{T}, 2}, phase_xys::Array{Complex{T}, 3}, k_x::Vector{T}, k_y::Vector{T}) where{T}

    for n in 1:size(poses, 1)
        x, y, z = poses[n]
        revise_phase_neg!((@view phase_xs[:, n]), (@view phase_ys[:, n]), k_x, k_y, x - L[1]/2, y - L[2]/2)
    end

    for n in 1:size(phase_xys, 3), j in 1:size(phase_xys, 2), i in 1:size(phase_xys, 1)
        phase_xys[i, j, n] = phase_xs[i, n] * phase_ys[j, n] * qs[n]
    end

    return nothing
end

# precompute the inverse matrix
function inverse_mat(N_grid::NTuple{3, Int}, L_z::T, k_x::Vector{T}, k_y::Vector{T}) where{T}
    N_x, N_y, N_z = N_grid
    inverse_mat = zeros(T, N_z + 2, N_z + 2, 2 * N_x + 1, 2 * N_y + 1)

    scale = L_z / T(2)

    A = zeros(T, N_z, N_z + 2)
    A[2,2] = 1 / 8 * scale^2
    A[2,4] = - 1 / 8 * scale^2
    for s in 3:N_z
        k = s - 1
        A[s, s - 2] = scale^2 / (2 * k * 2(k - 1))
        A[s, s] = - (1 / 2(k - 1) + 1 / 2(k + 1)) * scale^2 / 2k
        A[s, s + 2] = scale^2 / (2k * 2(k+1))
    end

    B = zeros(T, N_z, 2)
    B[1, 1] = T(2)
    B[2, 2] = scale

    C = zeros(T, 2, N_z + 2)
    for s = 3:N_z
        k = s - 1
        C[1, s-2] += scale^2 * (1/(2*k)) * (1/(2*(k-1))) * (-1)^k
        C[1, s] -= scale^2 * ((1/(2*k)) * (1/(2*(k-1))) + (1/(2*k)) * (1/(2*(k+1)))) * (-1)^k
        C[1, s+2] += scale^2 * (1/(2*k)) * (1/(2*(k+1))) * (-1)^k
        
        C[2, s-2] += scale^2 * (1/(2*k)) * (1/(2*(k-1)))
        C[2, s] -= scale^2 * ((1/(2*k)) * (1/(2*(k-1))) + (1/(2*k)) * (1/(2*(k+1))))
        C[2, s+2] += scale^2 * (1/(2*k)) * (1/(2*(k+1)))
    end
    C[1, 2] += 1/8 * scale^2 * (-1)
    C[1, 4] -= 1/8 * scale^2 * (-1)
    C[2, 2] += 1/8 * scale^2
    C[2, 4] -= 1/8 * scale^2

    D = zeros(T, 2, 2)
    D[1, 1] = T(1)
    D[1, 2] = - scale
    D[2, 1] = T(1)
    D[2, 2] = scale

    divisor = zeros(T, N_z + 2, N_z + 2)
    for i in 1:2N_x + 1, j in 1:2N_y + 1
        k = sqrt(k_x[i]^2 + k_y[j]^2)
        if k != zero(T)
            nu = - k^2
            divisor[1:N_z, 1:N_z] = A[1:N_z, 1:N_z] .* nu .+ T.(I(N_z))
            divisor[1:N_z, N_z+1:N_z+2] = nu .* B[1:N_z, 1:2]
            divisor[N_z+1:N_z+2, 1:N_z] = C[1:2, 1:N_z]
            divisor[N_z+1:N_z+2, N_z+1:N_z+2] = D

            inverse_mat[:,:,i,j] = inv(divisor)
        end
    end

    return inverse_mat
end
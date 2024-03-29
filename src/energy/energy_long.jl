# N_x, N_y are the uniform grid numbers in the x and y directions
# N_z is the Chebyshev grid number in the z direction
# n_atoms is the number of particles
# M_num is the number of Gaussians M_mid+1:M

# This part will be done in the following steps:
# 1. Generate the inverse matrix, k matrix
# 2. interpolate the particles onto the grid (Fourier grids in xy and chebgrid in z)
# 3. Convert the grid in z into the Chebyshev series
# 4. solve the equation for the boundary condition
# 5. gather the energy

function energy_long_loop_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, N_grid::NTuple{3, Int}, M_mid::Int, 
    k_x::Vector{T}, k_y::Vector{T}, r_z::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}},
    us_mat::Array{Complex{T}, 3}, b_l::Array{Complex{T}, 2}, b_u::Array{Complex{T}, 2}, 
    rhs::Vector{Complex{T}}, sol::Vector{Complex{T}}, ivsm::Array{T, 4},
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, H_s::Array{Complex{T}, 3}, 
    uspara::USeriesPara{T}
    ) where{T}

    @assert M_mid ≤ length(uspara.sw)

    b_l, b_u = boundaries!(qs, poses, L, b_l, b_u, k_x, k_y, us_mat, phase_x, phase_y, L[3], uspara, M_mid)
    H_r = interpolate_nu_loop!(H_r, qs, poses, N_grid, L, k_x, k_y, phase_x, phase_y, r_z, us_mat, uspara, M_mid)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    H_s = solve_eqs!(rhs, sol, H_c, H_s, b_l, b_u, ivsm, L[3])
    E_k = gather_nu(qs, poses, L, k_x, k_y, phase_x, phase_y, H_s)

    @debug "long range energy, einsum, FFCT" E_k

    return E_k
end

function energy_long_einsum_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, M_mid::Int, 
    k_x::Vector{T}, k_y::Vector{T}, k_mat::Array{T, 2}, r_z::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}},
    phase_xs::Array{Complex{T}, 2}, phase_ys::Array{Complex{T}, 2}, phase_xys::Array{Complex{T}, 3},
    temp_ijlk::Array{Complex{T}, 4}, temp_ijl::Array{Complex{T}, 3}, size_dict::Dict{Char, Int64}, z_coef::Array{T, 3}, exp_coef::Array{T, 3}, 
    us_mat::Array{Complex{T}, 3}, b_l::Array{Complex{T}, 2}, b_u::Array{Complex{T}, 2}, 
    rhs::Vector{Complex{T}}, sol::Vector{Complex{T}}, ivsm::Array{T, 4},
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, H_s::Array{Complex{T}, 3}, 
    uspara::USeriesPara{T}
    ) where{T}

    @assert M_mid ≤ length(uspara.sw)

    b_l, b_u = boundaries!(qs, poses, L, b_l, b_u, k_x, k_y, us_mat, phase_x, phase_y, L[3], uspara, M_mid)
    H_r = interpolate_nu_einsum!(H_r, qs, poses, L, k_x, k_y, k_mat, phase_xs, phase_ys, phase_xys, z_coef, exp_coef, temp_ijlk, temp_ijl, r_z, us_mat, uspara, M_mid, size_dict)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    H_s = solve_eqs!(rhs, sol, H_c, H_s, b_l, b_u, ivsm, L[3])
    E_k = gather_nu(qs, poses, L, k_x, k_y, phase_x, phase_y, H_s)

    @debug "long range energy, loop, FFCT" E_k

    return E_k
end

function energy_long_cheb_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, 
    M_mid::Int, 
    k_x::Vector{T}, k_y::Vector{T}, r_z::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}},
    us_mat::Array{Complex{T}, 3}, b_l::Array{Complex{T}, 2}, b_u::Array{Complex{T}, 2}, 
    rhs::Vector{Complex{T}}, sol::Vector{Complex{T}}, ivsm::Array{T, 4},
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, H_s::Array{Complex{T}, 3}, 
    uspara::USeriesPara{T}, cheb_mat::Array{ChebPoly{1, T, T}, 2}
    ) where{T}

    @assert M_mid ≤ length(uspara.sw)

    b_l, b_u = boundaries!(qs, poses, L, b_l, b_u, k_x, k_y, us_mat, phase_x, phase_y, L[3], uspara, M_mid)
    H_r = interpolate_nu_cheb!(H_r,qs, poses, L, k_x, k_y, phase_x, phase_y, r_z, cheb_mat)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    H_s = solve_eqs!(rhs, sol, H_c, H_s, b_l, b_u, ivsm, L[3])
    E_k = gather_nu(qs, poses, L, k_x, k_y, phase_x, phase_y, H_s)

    @debug "long range energy, cheb, FFCT" E_k

    return E_k
end

function energy_long_0(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, M_mid::Int, 
    z::Vector{T}, sort_z::Vector{Int}, 
    uspara::USeriesPara{T}, soepara::SoePara{Complex{T}}
    ) where{T}

    @assert M_mid ≤ length(uspara.sw)

    revise_z!(z, sort_z, poses)
    E_0 = zeroth_order(qs, z, soepara, uspara, sort_z, L, M_mid)

    return E_0
end

function energy_long(interaction::FSSoGInteraction{T}) where{T}

    E_k = energy_long_cheb_k(interaction.charge, interaction.position, interaction.L, interaction.M_mid, interaction.k_x, interaction.k_y, interaction.r_z, interaction.phase_x, interaction.phase_y, interaction.us_mat, interaction.b_l, interaction.b_u, interaction.rhs, interaction.sol, interaction.ivsm, interaction.H_r, interaction.H_c, interaction.H_s, interaction.uspara, interaction.cheb_mat)
    E_0 = energy_long_0(interaction.charge, interaction.position, interaction.L, interaction.M_mid, interaction.z, interaction.sort_z, interaction.uspara, interaction.soepara)

    return (E_k + E_0) / (4π * interaction.ϵ)
end

function energy_long_Q2D_direct_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, 
    M_mid::Int, 
    k_x::Vector{T}, k_y::Vector{T}, r_z::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}},
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3},
    uspara::USeriesPara{T}) where{T}

    @assert M_mid ≤ length(uspara.sw)

    H_r = interpolate_Q2D_direct!(H_r, qs, poses, L, k_x, k_y, phase_x, phase_y, r_z, uspara, M_mid)
    H_c = real2Cheb_Q2D!(H_r, H_c, r_z, L[3])
    E_k = gather_Q2D(qs, poses, L, k_x, k_y, phase_x, phase_y, H_c)

    @debug "long range energy, FFCT, direct interpolate" E_k

    return E_k
end
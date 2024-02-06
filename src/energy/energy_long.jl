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

function energy_long(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, M_mid::Int, 
    k_x::Vector{T}, k_y::Vector{T}, r_z::Vector{T}, 
    z::Vector{T}, sort_z::Vector{Int}, 
    us_mat::Array{T, 3}, b_l::Array{Complex{T}, 2}, b_u::Array{Complex{T}, 2}, 
    rhs::Vector{Complex{T}}, sol::Vector{Complex{T}}, ivsm::Array{T, 4},
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, H_s::Array{Complex{T}, 3}, 
    uspara::USeriesPara{T}, soepara::SoePara{Complex{T}}
    ) where{T}

    @assert M_mid ≤ length(uspara.sw)

    b_l, b_u = boundaries!(qs, poses, b_l, b_u, k_x, k_y, L[3], uspara, M_mid)
    H_r = interpolate_nu!(qs, poses, k_x, k_y, r_z, us_mat, H_r, uspara, M_mid)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    H_s = solve_eqs!(rhs, sol, H_c, H_s, b_l, b_u, ivsm, L[3])
    E_k = gather_nu(qs, poses, L, k_x, k_y, H_s)

    revise_z!(z, sort_z, poses)
    E_0 = zeroth_order(qs, z, soepara, uspara, sort_z, L, M_mid)

    @debug "long range energy, FFCT" Ek, E0

    return E_k + E_0
end
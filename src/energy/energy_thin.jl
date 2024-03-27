function energy_thin_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}},
    qs_new::Vector{T}, poses_new::Vector{NTuple{2, T}},
    L::NTuple{3, T}, r_z::Vector{T}, 
    b_l::Array{Complex{T}, 2}, b_u::Array{Complex{T}, 2}, 
    rhs::Vector{Complex{T}}, sol::Vector{Complex{T}}, ivsm::Array{T, 4},
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, H_s::Array{Complex{T}, 3}, 
    uspara::USeriesPara{T}, 
    gridinfo::GridInfo{2, T}, gridboxs::Array{GridBox{2, T}, 2}, cheb_coefs::NTuple{2, ChebCoef{T}}, scalefactors::Vector{ScalingFactor{T}}
    ) where{T}

    @assert M_mid ≤ length(uspara.sw)

    b_l, b_u = boundaries_thin!()
    gridboxs = interpolate_thin!(qs, poses, qs_new, poses_new, L, gridinfo, gridboxs, cheb_coefs, r_z)
    gridboxs = thin_scale!(gridboxs, scalefactors)
    H_r = sum_gridboxs!(gridboxs, H_r)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    H_c = ifft2d!(H_c)
    H_s = solve_eqs!(rhs, sol, H_c, H_s, b_l, b_u, ivsm, L[3])
    E_k = gather_thin(qs, poses, L, H_s, gridinfo, cheb_value, cheb_coefs)

    @debug "long range energy, thin" E_k

    return E_k
end
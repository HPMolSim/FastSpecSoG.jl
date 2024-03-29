function energy_thin_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}},
    qs_new::Vector{T}, poses_new::Vector{NTuple{2, T}},
    L::NTuple{3, T}, r_z::Vector{T}, 
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3},
    gridinfo::GridInfo{2, T}, gridboxs::Array{GridBox{2, T}, 2}, cheb_coefs::NTuple{2, ChebCoef{T}}, scalefactors::Vector{ScalingFactor{2, T}}, cheb_value::Vector{Array{T, 1}}
    ) where{T}

    gridboxs = interpolate_thin!(qs, poses, qs_new, poses_new, L, gridinfo, gridboxs, cheb_coefs, r_z)
    H_r = scale_thin!(gridboxs, scalefactors, H_r)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    H_c = ifft2d!(H_c)
    E_k = gather_thin(qs, poses, L, H_c, gridinfo, cheb_value, cheb_coefs)

    @debug "long range energy, thin" E_k

    return E_k
end
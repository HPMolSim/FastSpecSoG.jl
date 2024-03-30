function energy_long_thin_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}},
    L::NTuple{3, T}, r_z::Vector{T}, 
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3},
    gridinfo::GridInfo{2, T}, pad_grids::Vector{Array{Complex{T}, 3}}, scalefactors::Vector{ScalingFactor{2, T}}, 
    cheb_coefs::NTuple{2, ChebCoef{T}}, cheb_value::Vector{Array{T, 1}}
    ) where{T}

    gridboxs = interpolate_thin!(qs, poses, pad_grids, gridinfo, cheb_value, cheb_coefs, r_z, L[3])
    piecewise_fft!.(pad_grids)
    H_r = scale_thin!(pad_grids, scalefactors, H_r)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    H_c = piecewise_ifft!(H_c)
    E_k = gather_thin(qs, poses, L, H_c, gridinfo, cheb_value, cheb_coefs)

    @debug "long range energy, thin" E_k

    return E_k
end
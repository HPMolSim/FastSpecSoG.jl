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

    E_k = energy_long_cheb_k(interaction.charge, interaction.position, interaction.L, interaction.k_x, interaction.k_y, interaction.r_z, interaction.phase_x, interaction.phase_y, interaction.H_r, interaction.H_c, interaction.cheb_mat)

    E_0 = energy_long_0(interaction.charge, interaction.position, interaction.L, interaction.M_mid, interaction.z, interaction.sort_z, interaction.uspara, interaction.soepara)

    return (E_k + E_0) / (4π * interaction.ϵ)
end

function energy_long(interaction::FSSoGThinInteraction{T}) where{T}

    E_k = energy_long_thin_k(interaction.charge, interaction.position, interaction.L, interaction.r_z, interaction.H_r, interaction.H_c, interaction.gridinfo, interaction.pad_grids, interaction.scalefactors, interaction.cheb_coefs, interaction.cheb_value)

    E_0 = energy_long_0(interaction.charge, interaction.position, interaction.L, 0, interaction.z, interaction.sort_z, interaction.uspara, interaction.soepara)
    return (E_k + E_0) / (4π * interaction.ϵ)
end

function energy_long_loop_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, 
    M_mid::Int, 
    k_x::Vector{T}, k_y::Vector{T}, r_z::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}},
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3},
    uspara::USeriesPara{T}) where{T}

    @assert M_mid ≤ length(uspara.sw)

    H_r = interpolate_nu_loop!(H_r, qs, poses, L, k_x, k_y, phase_x, phase_y, r_z, uspara, M_mid)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    E_k = gather_nu(qs, poses, L, k_x, k_y, phase_x, phase_y, H_c)

    @debug "long range energy, cube, loop" E_k

    return E_k
end

function energy_long_cheb_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, 
    k_x::Vector{T}, k_y::Vector{T}, r_z::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}},
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, cheb_mat::Array{ChebPoly{1, T, T}, 2}) where{T}

    H_r = interpolate_nu_cheb!(H_r, qs, poses, L, k_x, k_y, phase_x, phase_y, r_z, cheb_mat)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    E_k = gather_nu(qs, poses, L, k_x, k_y, phase_x, phase_y, H_c)

    @debug "long range energy, cube, cheb" E_k

    return E_k
end

function energy_long_thin_k(
    qs::Vector{T}, poses::Vector{NTuple{3, T}},
    L::NTuple{3, T}, r_z::Vector{T}, 
    H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3},
    gridinfo::GridInfo{2, T}, pad_grids::Vector{Array{Complex{T}, 3}}, scalefactors::Vector{ScalingFactor{2, T}}, 
    cheb_coefs::NTuple{2, ChebCoef{T}}, cheb_value::Vector{Array{T, 1}}
    ) where{T}

    interpolate_thin!(qs, poses, pad_grids, gridinfo, cheb_value, cheb_coefs, r_z, L[3])
    piecewise_fft!.(pad_grids)
    H_r = scale_thin!(pad_grids, scalefactors, H_r)
    H_c = real2Cheb!(H_r, H_c, r_z, L[3])
    H_c = piecewise_ifft!(H_c)
    E_k = gather_thin(qs, poses, L, H_c, gridinfo, cheb_value, cheb_coefs)

    @debug "long range energy, thin" E_k

    return E_k
end
# precomputings

@inline @inbounds function set_zeros!(a::AbstractArray{T}) where{T <: Number}
    return fill!(a, zero(T))
end

# precompute the parameters to be used
function long_paras_gen(L::NTuple{3, T}, N_grid::NTuple{3, Int}, n_atoms::Int) where{T}
    L_x, L_y, L_z = L
    N_x, N_y, N_z = N_grid
    k_x = Vector{T}([2π * i / L_x for i in -N_x:N_x])
    k_y = Vector{T}([2π * i / L_y for i in -N_y:N_y])
    r_z = Vector{T}([L_z / 2 * ( 1.0 - cos((2i - 1)*π / 2N_z) ) for i in 1:N_z])

    # xy k space grid, z real space grid
    H_r = zeros(Complex{T}, 2N_x + 1, 2N_y + 1, N_z)
    # xy k space grid, z ChebCoefs
    H_c = zeros(Complex{T}, 2N_x + 1, 2N_y + 1, N_z)

    phase_x = zeros(Complex{T}, 2N_x + 1)
    phase_y = zeros(Complex{T}, 2N_y + 1)

    sort_z = zeros(Int, n_atoms)
    z = zeros(T, n_atoms)
    
    return k_x, k_y, r_z, H_r, H_c, phase_x, phase_y,  sort_z, z
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

function thin_grids_gen(N_real::NTuple{2, Int}, R_z::Int, w::NTuple{2, Int}, β::NTuple{2, T}, L::NTuple{3, T}, cheb_order::Int, uspara::USeriesPara{T}, Taylor_Q::Int) where{T}
    periodicity = (true, true)
    extra_pad = (0, 0)

    gridinfo = GridInfo(N_real, w, periodicity, extra_pad, (L[1], L[2]))
    pad_grids = [zeros(Complex{T}, N_real[1], N_real[2], R_z) for i in 1:Taylor_Q]

    f_window = [x -> Wkb(x, (w[i] + 0.5) * gridinfo.h[i], β[i]) for i in 1:2]
    cheb_coefs = tuple([ChebCoef(f_window[i], gridinfo.h[i], w[i], cheb_order) for i in 1:2]...)

    F_f_window = [x -> FWkb(x, (w[i] + 0.5) * gridinfo.h[i], β[i]) for i in 1:2]

    k_x = gridinfo.k[1]
    k_y = gridinfo.k[2]
    taylor_mats = TaylorUSeries(k_x, k_y, L[3], uspara, 0, Taylor_Q)

    func_scale = (k_x, k_y) -> (F_f_window[1](k_x) * F_f_window[2](k_y))^(-2)
    sf0 = ScalingFactor(func_scale, gridinfo)
    scalefactors = [ScalingFactor(func_scale, sf0.factors .* taylor_mats[i]) for i in 1:Taylor_Q]

    return (gridinfo, pad_grids, cheb_coefs, scalefactors)
end

function thin_paras_gen(N_real::NTuple{2, Int}, R_z::Int, w::NTuple{2, Int}, β::NTuple{2, T}, L::NTuple{3, T}, cheb_order::Int, uspara::USeriesPara{T}, Taylor_Q::Int) where{T}
    gridinfo, pad_grids, cheb_coefs, scalefactors = thin_grids_gen(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)

    H_r = zeros(Complex{T}, N_real[1], N_real[2], R_z)
    H_c = zeros(Complex{T}, N_real[1], N_real[2], R_z)

    cheb_value = [zeros(T, size(cheb_coefs[i].coef, 1) - 1) for i in 1:2]
    r_z = Vector{T}([L[3] / 2 * ( 1.0 - cos((2i - 1)*π / 2R_z) ) for i in 1:R_z])

    return gridinfo, pad_grids, cheb_coefs, scalefactors, H_r, H_c, cheb_value, r_z
end
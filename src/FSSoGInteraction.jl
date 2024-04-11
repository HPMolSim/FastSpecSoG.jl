function k_set_2D(k_c::T, L::NTuple{3,T}) where{T}
    L_x, L_y, L_z = L

    mx_max = ceil(Int, k_c * L_x / 2π) + 1
    my_max = ceil(Int, k_c * L_y / 2π) + 1

    k_set = Vector{NTuple{3, T}}()

    for m_x in - mx_max : mx_max
        for m_y in - my_max : my_max
            k_x = m_x * 2π / L_x
            k_y = m_y * 2π / L_y
            k = sqrt(k_x^2 + k_y^2)
            if k <= k_c
                push!(k_set, (k_x, k_y, k))
            end
        end
    end

    return k_set
end


function FSSoG_naive(L::NTuple{3, T}, n_atoms::Int, r_c::T, k_c::T, b::T, σ::T, ω::T, M::Int; ϵ::T = 1.0) where{T}
    uspara = USeriesPara(b, σ, ω, M)
    k_set = k_set_2D(k_c, L)
    return FSSoG_naive{T}(b, σ, ω, M, r_c, k_c, ϵ, L, k_set, uspara, n_atoms)
end

function FSSoG_naive(L::NTuple{3, T}, n_atoms::Int, r_c::T, k_c::T; preset::Int = 1, ϵ::T = 1.0) where{T}

    @assert preset ≤ length(preset_parameters)
    
    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])
    uspara = USeriesPara(b, σ, ω, M)
    k_set = k_set_2D(k_c, L)
    return FSSoG_naive{T}(b, σ, ω, M, r_c, k_c, ϵ, L, k_set, uspara, n_atoms)
end

function FSSoGInteraction(
    b::T, σ::T, ω::T, M::Int,
    L::NTuple{3, T}, n_atoms::Int,
    r_c::T, Q_short::Int, r_min::T,
    N_real_mid::NTuple{3, Int}, w::NTuple{3, Int}, β::NTuple{3, T}, extra_pad_ratio::Int, cheb_order::Int, M_mid::Int, 
    N_grid_long::NTuple{3, Int}, Q_long::Int, 
    Rz_0::Int, Q_0::Int; 
    ϵ::T = one(T)) where{T}

    uspara = USeriesPara(b, σ, ω, M)

    position = [tuple(zeros(T, 3)...) for i in 1:n_atoms]
    charge = zeros(T, n_atoms)

    boundary = Q2dBoundary(L...)

    uspara_cheb = Es_USeries_Cheb(uspara, r_min, r_c, Q_short)
    F0 = F0_cal(b, σ, ω, M)

    gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real_mid, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
    
    # some of these term are for other methods, not used
    k_x, k_y, r_z, H_r, H_c, phase_x, phase_y = long_paras_gen(L, N_grid_long)
    cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, Q_long)

    r_z0, grids0, chebcoefs0 = zero_paras_gen(L[3], Rz_0)
    chebuseries = ChebUSeries_0(L[3], uspara, M_mid, Q_0)

    return FSSoGInteraction{T}(b, σ, ω, M, ϵ, L, boundary, n_atoms, uspara, position, charge, uspara_cheb, r_c, F0, gridinfo, gridbox, cheb_coefs, scalefactor, M_mid, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, cheb_mat, Q_0, r_z0, chebcoefs0, grids0, chebuseries)
end

function FSSoGInteraction(
    L::NTuple{3, T}, n_atoms::Int,
    r_c::T, Q_short::Int, r_min::T,
    N_real_mid::NTuple{3, Int}, w::NTuple{3, Int}, β::NTuple{3, T}, extra_pad_ratio::Int, cheb_order::Int, M_mid::Int, 
    N_grid_long::NTuple{3, Int}, Q_long::Int, 
    Rz_0::Int, Q_0::Int; 
    preset::Int = 1, ϵ::T = one(T)) where{T}

    @assert preset ≤ length(preset_parameters)
    
    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])
    
    return FSSoGInteraction(b, σ, ω, M, L, n_atoms, r_c, Q_short, r_min, N_real_mid, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid_long, Q_long, Rz_0, Q_0, ϵ = ϵ)
end

function FSSoGThinInteraction(
    b::T, σ::T, ω::T, M::Int,
    L::NTuple{3, T}, n_atoms::Int,
    r_c::T, Q_short::Int, r_min::T,
    N_real::NTuple{2, Int}, R_z::Int, w::NTuple{2, Int}, β::NTuple{2, T}, cheb_order::Int, Taylor_Q::Int,
    Rz_0::Int, Q_0::Int; 
    ϵ::T = one(T)) where{T}

    uspara = USeriesPara(b, σ, ω, M)
    position = [tuple(zeros(T, 3)...) for i in 1:n_atoms]
    charge = zeros(T, n_atoms)

    boundary = Q2dBoundary(L...)

    uspara_cheb = Es_USeries_Cheb(uspara, r_min, r_c, Q_short)
    F0 = F0_cal(b, σ, ω, M)

    gridinfo, pad_grids, cheb_coefs, scalefactors, H_r, H_c, cheb_value, r_z = thin_paras_gen(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)

    r_z0, grids0, chebcoefs0 = zero_paras_gen(L[3], Rz_0)
    chebuseries = ChebUSeries_0(L[3], uspara, 0, Q_0)

    return FSSoGThinInteraction{T}(b, σ, ω, M, ϵ, L, boundary, n_atoms, uspara, position, charge, uspara_cheb, r_c, F0, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value, Q_0, r_z0, chebcoefs0, grids0, chebuseries)
end

function FSSoGThinInteraction(
    L::NTuple{3, T}, n_atoms::Int,
    r_c::T, Q_short::Int, r_min::T,
    N_real::NTuple{2, Int}, R_z::Int, w::NTuple{2, Int}, β::NTuple{2, T}, cheb_order::Int, Taylor_Q::Int,
    Rz_0::Int, Q_0::Int; 
    preset::Int = 1, ϵ::T = one(T)) where{T}

    @assert preset ≤ length(preset_parameters)
    
    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])
    
    return FSSoGThinInteraction(b, σ, ω, M, L, n_atoms, r_c, Q_short, r_min, N_real, R_z, w, β, cheb_order, Taylor_Q, Rz_0, Q_0, ϵ = ϵ)
end
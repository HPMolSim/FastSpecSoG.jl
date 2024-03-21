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
    N_grid_long::NTuple{3, Int}, Q_long::Int, soepara::SoePara{Complex{T}}; 
    ϵ::T = one(T)) where{T}

    uspara = USeriesPara(b, σ, ω, M)

    position = [tuple(zeros(T, 3)...) for i in 1:n_atoms]
    charge = zeros(T, n_atoms)

    boundary = Q2dBoundary(L...)

    uspara_cheb = Es_USeries_Cheb(uspara, r_min, r_c, Q_short)
    F0 = F0_cal(b, σ, ω, M)

    gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real_mid, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)
    
    # some of these term are for other methods, not used
    k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid_long, uspara, M_mid, n_atoms)
    cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, Q_long)

    return FSSoGInteraction{T}(b, σ, ω, M, ϵ, L, boundary, n_atoms, uspara, position, charge, uspara_cheb, r_c, F0, gridinfo, gridbox, cheb_coefs, scalefactor, M_mid, k_x, k_y, r_z, phase_x, phase_y, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, cheb_mat, z, sort_z, soepara)
end

function FSSoGInteraction(
    L::NTuple{3, T}, n_atoms::Int,
    r_c::T, Q_short::Int, r_min::T,
    N_real_mid::NTuple{3, Int}, w::NTuple{3, Int}, β::NTuple{3, T}, extra_pad_ratio::Int, cheb_order::Int, M_mid::Int, 
    N_grid_long::NTuple{3, Int}, Q_long::Int, soepara::SoePara{Complex{T}}; 
    preset::Int = 1, ϵ::T = one(T)) where{T}

    @assert preset ≤ length(preset_parameters)
    
    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])
    
    return FSSoGInteraction(b, σ, ω, M, L, n_atoms, r_c, Q_short, r_min, N_real_mid, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid_long, Q_long, soepara, ϵ = ϵ)
end
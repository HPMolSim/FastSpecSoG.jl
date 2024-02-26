@inline function Fourier_k(i::Int, N::Int, L::T) where{T}
    return T(2π * (i - 1 - N) / L)
end


@inbounds function revise_z_coef!(z_coef::Array{T, 3}, exp_coef::Array{T, 3}, r_z::Vector{T}, poses::Vector{NTuple{3, T}}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    for n in 1:size(poses, 1)
        x, y, z = poses[n]
        for k in 1:size(r_z, 1)
            r_zk = r_z[k]
            for l in M_mid + 1:length(uspara.sw)
                sl, wl = uspara.sw[l]
                temp_z = (z - r_zk)^2 / sl^2
                temp_exp = exp(-temp_z)
                exp_coef[l - M_mid, k, n] = temp_exp
                z_coef[l - M_mid, k, n] = (T(2) - T(4) * temp_z) * temp_exp * sl^2
            end
        end
    end

    return nothing
end


"""
H1[i, j, k] := phase_xys[i, j, n] * z_coef[l, k, n] * us_mat[i, j, l]
H2[i, j, k] := phase_xys[i, j, n] * exp_coef[l, k, n] * k_mat[i, j] * us_mat[i, j, l]
"""

@inbounds function interpolate_nu_einsum!(
    H_r::Array{Complex{T}, 3},
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, 
    k_x::Vector{T}, k_y::Vector{T}, k_mat::Array{T, 2}, 
    phase_xs::Array{Complex{T}, 2}, phase_ys::Array{Complex{T}, 2}, phase_xys::Array{Complex{T}, 3}, 
    z_coef::Array{T, 3}, exp_coef::Array{T, 3}, temp_ijlk::Array{Complex{T}, 4}, temp_ijl::Array{Complex{T}, 3},
    r_z::Vector{T}, us_mat::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int,  
    size_dict::Dict{Char, Int64}) where{T}

    revise_phase_neg_all!(qs, poses, phase_xs, phase_ys, phase_xys, k_x, k_y)
    revise_z_coef!(z_coef, exp_coef, r_z, poses, uspara, M_mid)

    # H1
    einsum!(ein"ijn, lkn -> ijlk", (phase_xys, z_coef), temp_ijlk, true, false, size_dict)
    einsum!(ein"ijlk, ijl -> ijk", (temp_ijlk, us_mat), H_r, true, false, size_dict)

    # H2
    einsum!(ein"ijn, lkn -> ijlk", (phase_xys, exp_coef), temp_ijlk, true, false, size_dict)
    einsum!(ein"ij, ijl -> ijl", (k_mat, us_mat), temp_ijl, true, false, size_dict)
    einsum!(ein"ijlk, ijl -> ijk", (temp_ijlk, temp_ijl), H_r, true, true, size_dict)

    H_r .*= π

    return H_r
end

@inbounds function interpolate_nu_einsum_non_inplace!(
    H_r::Array{Complex{T}, 3},
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, 
    k_x::Vector{T}, k_y::Vector{T}, k_mat::Array{T, 2}, 
    phase_xs::Array{Complex{T}, 2}, phase_ys::Array{Complex{T}, 2}, phase_xys::Array{Complex{T}, 3}, 
    z_coef::Array{T, 3}, exp_coef::Array{T, 3}, 
    r_z::Vector{T}, us_mat::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int,  
    size_dict::Dict{Char, Int64}) where{T}

    revise_phase_neg_all!(qs, poses, phase_xs, phase_ys, phase_xys, k_x, k_y)
    revise_z_coef!(z_coef, exp_coef, r_z, poses, uspara, M_mid)

    temp_ijlk = ein"ijn, lkn -> ijlk"(phase_xys, z_coef)
    H_1 = ein"ijlk, ijl -> ijk"(temp_ijlk, us_mat)

    temp_ijlk = ein"ijn, lkn -> ijlk"(phase_xys, exp_coef)
    temp_ijl = ein"ij, ijl -> ijl"(k_mat, us_mat)
    H_2 = ein"ijlk, ijl -> ijk"(temp_ijlk, temp_ijl)

    H_r = π .* (H_1 + H_2)

    return H_r
end

@inbounds function interpolate_nu_loop_single!(
    H_r::Array{Complex{T}, 3},
    q::T, pos::NTuple{3, T}, 
    N::NTuple{3, Int}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, 
    r_z::Vector{T}, us_mat::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    x, y, z = pos
    # revise_phase_neg!(phase_x, phase_y, k_x, k_y, x, y)

    for k in 1:size(H_r, 3)
        r_zk = r_z[k]
        for l in M_mid + 1:length(uspara.sw)
            sl, wl = uspara.sw[l]

            for j in 1:size(H_r, 2)
                k_yj = k_y[j]
                for i in 1:size(H_r, 1)
                    k_xi = k_x[i]
                    phase = exp(-T(1)im * (k_xi * (x - L[1] / 2) + k_yj * (y - L[2] / 2)))
                    H_r[i, j, k] += q * π * wl * phase * (T(2) - T(4) * (z - r_zk)^2 / sl^2 + (k_xi^2 + k_yj^2) * sl^2) * exp(- (z - r_zk)^2 / sl^2) * exp(- sl^2 * (k_xi^2 + k_yj^2) / 4)
                end
            end
        end
    end

    return H_r
end

function interpolate_nu_loop!(
    H_r::Array{Complex{T}, 3},
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, 
    N::NTuple{3, Int}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, 
    r_z::Vector{T}, us_mat::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    set_zeros!(H_r)
    for i in 1:length(qs)
        interpolate_nu_loop_single!(H_r, qs[i], poses[i], N, L, k_x, k_y, phase_x, phase_y, r_z, us_mat, uspara, M_mid)
    end

    return H_r
end

@inbounds function real2Cheb!(H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

    set_zeros!(H_c)
    N_z = size(H_r, 3)

    for k in 1:N_z, l in 1:N_z
        cheb_temp = chebpoly(k - 1, r_z[l] - L_z / T(2), L_z / T(2))
        for j in 1:size(H_r, 2), i in 1:size(H_r, 1)
            H_c[i, j, k] += 2 / N_z * H_r[i, j, l] * cheb_temp
        end
    end

    return H_c
end
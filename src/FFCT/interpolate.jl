@inbounds function revise_z_coef!(z_coef::Array{T, 3}, exp_coef::Array{T, 3}, r_z::Vector{T}, poses::Vector{NTuple{3, T}}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    for n in 1:size(poses, 1)
        x, y, z = poses[n]
        for k in 1:size(r_z, 1)
            r_zk = r_z[k]
            for l in M_mid + 1:length(uspara.sw)
                sl, wl = uspara.sw[l]
                temp_z = (z - r_zk)^2 / sl^2
                temp_exp = wl * exp(-temp_z)
                exp_coef[l - M_mid, k, n] = temp_exp * sl^2
                z_coef[l - M_mid, k, n] = (T(2) - T(4) * temp_z) * temp_exp
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
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T},
    k_x::Vector{T}, k_y::Vector{T}, k_mat::Array{T, 2}, 
    phase_xs::Array{Complex{T}, 2}, phase_ys::Array{Complex{T}, 2}, phase_xys::Array{Complex{T}, 3}, 
    z_coef::Array{T, 3}, exp_coef::Array{T, 3}, temp_ijlk::Array{Complex{T}, 4}, temp_ijl::Array{Complex{T}, 3},
    r_z::Vector{T}, us_mat::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int, size_dict::Dict{Char, Int64}) where{T}

    revise_phase_neg_all!(qs, poses, L, phase_xs, phase_ys, phase_xys, k_x, k_y)
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
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T},
    k_x::Vector{T}, k_y::Vector{T}, k_mat::Array{T, 2}, 
    phase_xs::Array{Complex{T}, 2}, phase_ys::Array{Complex{T}, 2}, phase_xys::Array{Complex{T}, 3}, 
    z_coef::Array{T, 3}, exp_coef::Array{T, 3}, 
    r_z::Vector{T}, us_mat::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    revise_phase_neg_all!(qs, poses, L, phase_xs, phase_ys, phase_xys, k_x, k_y)
    revise_z_coef!(z_coef, exp_coef, r_z, poses, uspara, M_mid)

    @ein H1[i, j, k] := phase_xys[i, j, n] * z_coef[l, k, n] * us_mat[i, j, l]
    @ein H2[i, j, k] := phase_xys[i, j, n] * exp_coef[l, k, n] * k_mat[i, j] * us_mat[i, j, l]

    H_r = π .* (H1 .+ H2)

    return H_r
end

@inbounds function interpolate_nu_cheb_single!(
    H_r::Array{Complex{T}, 3},
    q::T, pos::NTuple{3, T}, L::NTuple{3, T},
    k_x::Vector{T}, k_y::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, 
    r_z::Vector{T}, cheb_mat::Array{ChebPoly{1, T, T}, 2}) where{T}

    x, y, z = pos
    revise_phase_neg!(phase_x, phase_y, k_x, k_y, x - L[1] / 2, y - L[2] / 2)

    
    for j in 1:size(H_r, 2)
        phase_yj = phase_y[j]
        for i in 1:size(H_r, 1)
            phase_xi = phase_x[i]
            phase = phase_xi * phase_yj
            f_cheb = cheb_mat[i, j]
            for k in 1:size(H_r, 3)
                r_zk = r_z[k]
                H_r[i, j, k] += q * phase * f_cheb(abs(z - r_zk))
            end
        end
    end

    return H_r
end

@inbounds function interpolate_nu_cheb!(
    H_r::Array{Complex{T}, 3},
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T},
    k_x::Vector{T}, k_y::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, 
    r_z::Vector{T}, cheb_mat::Array{ChebPoly{1, T, T}, 2}) where{T}

    set_zeros!(H_r)
    for i in 1:length(qs)
        interpolate_nu_cheb_single!(H_r, qs[i], poses[i], L, k_x, k_y, phase_x, phase_y, r_z, cheb_mat)
    end

    return H_r
end

@inbounds function interpolate_nu_loop_single!(
    H_r::Array{Complex{T}, 3},
    q::T, pos::NTuple{3, T}, 
    N::NTuple{3, Int}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, 
    r_z::Vector{T}, us_mat::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    x, y, z = pos
    revise_phase_neg!(phase_x, phase_y, k_x, k_y, x - L[1] / 2, y - L[2] / 2)

    for k in 1:size(H_r, 3)
        r_zk = r_z[k]
        for l in M_mid + 1:length(uspara.sw)
            sl, wl = uspara.sw[l]

            for j in 1:size(H_r, 2)
                k_yj = k_y[j]
                for i in 1:size(H_r, 1)
                    k_xi = k_x[i]
                    phase = phase_x[i] * phase_y[j]
                    us = us_mat[i, j, l - M_mid]
                    H_r[i, j, k] += q * π * wl * phase * (T(2) - T(4) * (z - r_zk)^2 / sl^2 + (k_xi^2 + k_yj^2) * sl^2) * exp(- (z - r_zk)^2 / sl^2) * us
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

@inbounds function real2Cheb_Q2D!(H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

    set_zeros!(H_c)
    N_z = size(H_r, 3)

    for k in 1:N_z, l in 1:N_z
        cheb_temp = k == 1 ? 0.5 : chebpoly(k - 1, r_z[l] - L_z / T(2), L_z / T(2))
        for j in 1:size(H_r, 2), i in 1:size(H_r, 1)
            H_c[i, j, k] += 2 / N_z * H_r[i, j, l] * cheb_temp
        end
    end

    return H_c
end

@inbounds function Cheb2real_Q2D!(H_c::Array{Complex{T}, 3}, H_r::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

    set_zeros!(H_r)
    N_z = size(H_c, 3)

    for k in 1:N_z, l in 1:N_z
        cheb_temp = chebpoly(k - 1, r_z[l] - L_z / T(2), L_z / T(2))
        for j in 1:size(H_c, 2), i in 1:size(H_c, 1)
            H_r[i, j, l] += H_c[i, j, k] * cheb_temp
        end
    end

    return H_r
end

function interpolate_Q2D_direct!(
    H_r::Array{Complex{T}, 3},
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T},
    k_x::Vector{T}, k_y::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, 
    r_z::Vector{T}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    set_zeros!(H_r)

    for n in 1:length(qs)
        x, y, z = poses[n]
        q = qs[n]
        revise_phase_neg!(phase_x, phase_y, k_x, k_y, x - L[1] / 2, y - L[2] / 2)

        for k in 1:size(H_r, 3)
            r_zk = r_z[k]
            for l in M_mid + 1:length(uspara.sw)
                sl, wl = uspara.sw[l]

                for j in 1:size(H_r, 2)
                    k_yj = k_y[j]
                    for i in 1:size(H_r, 1)
                        k_xi = k_x[i]
                        k2 = k_xi^2 + k_yj^2
                        phase = phase_x[i] * phase_y[j]
                        H_r[i, j, k] += q * π * wl * sl^2 * phase * exp(-sl^2 * k2 / 4) * exp(-(z - r_zk)^2 / sl^2)
                    end
                end
            end
        end
    end

    return H_r
end

function update_poses_new!(poses_new::Vector{NTuple{2, T}}, poses::Vector{NTuple{3, T}}) where{T}
    for i in 1:length(poses)
        poses_new[i] = (poses[i][1], poses[i][2])
    end
    return poses_new
end

function update_qs_new!(qs_new::Vector{T}, qs::Vector{T}, poses::Vector{NTuple{3, T}}, r_z::T, j::Int, Lz::T) where{T}
    for i in 1:length(qs)
        qs_new[i] = qs[i] * ((poses[i][3] - r_z) / Lz)^(2 * (j - 1))
    end
    return qs_new
end

function interpolate_thin!(
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, qs_new::Vector{T}, poses_new::Vector{NTuple{2, T}}, L::NTuple{3, T},
    gridinfo::GridInfo{2, T}, gridboxs::Array{GridBox{2, T}, 2}, cheb_coefs::NTuple{2, ChebCoef{T}},
    r_z::Vector{T}
    ) where{T}
    
    R_z, Taylor_Q = size(gridboxs)
    update_poses_new!(poses_new, poses)

    for i in 1:R_z
        r_zi = r_z[i]
        for j in 1:Taylor_Q
            update_qs_new!(qs_new, qs, poses, r_zi, j, L[3])
            interpolate!(qs_new, poses_new, gridinfo, gridboxs[i, j], cheb_coefs)
        end
    end

    return gridboxs
end
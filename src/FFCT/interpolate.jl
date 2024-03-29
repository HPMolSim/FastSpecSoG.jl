function interpolate_nu_loop!(
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
                        if !(k2 ≈ zero(T))
                            H_r[i, j, k] += q * π * wl * sl^2 * phase * exp(-sl^2 * k2 / 4) * exp(-(z - r_zk)^2 / sl^2)
                        end
                    end
                end
            end
        end
    end

    return H_r
end

function interpolate_nu_cheb!(
    H_r::Array{Complex{T}, 3},
    qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T},
    k_x::Vector{T}, k_y::Vector{T}, 
    phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, 
    r_z::Vector{T}, cheb_mat::Array{ChebPoly{1, T, T}, 2}) where{T}

    set_zeros!(H_r)

    for n in 1:length(qs)
        x, y, z = poses[n]
        q = qs[n]
        revise_phase_neg!(phase_x, phase_y, k_x, k_y, x - L[1] / 2, y - L[2] / 2)

        for k in 1:size(H_r, 3)
            r_zk = r_z[k]
            for j in 1:size(H_r, 2)
                k_yj = k_y[j]
                for i in 1:size(H_r, 1)
                    k_xi = k_x[i]
                    k2 = k_xi^2 + k_yj^2
                    if !(k2 ≈ zero(T))
                        phase = phase_x[i] * phase_y[j]
                        H_r[i, j, k] += q * phase * cheb_mat[i, j](abs(z - r_zk))
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

function ϕkz_direct(z0::T, qs::Vector{T}, poses::Vector{NTuple{3, T}}, k_x::T, k_y::T, uspara::USeriesPara{T}, M_mid::Int) where{T}

    s = zero(T)

    for n in 1:length(qs)
        x, y, z = poses[n]
        q = qs[n]

        for l in M_mid + 1:length(uspara.sw)
            sl, wl = uspara.sw[l]
            k2 = k_x^2 + k_y^2
            phase = exp( - T(1)im * (k_x * x + k_y * y))
            s += q * π * wl * sl^2 * phase * exp(-sl^2 * k2 / 4) * exp(-(z - z0)^2 / sl^2)
        end
    end

    return s
end
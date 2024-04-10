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

        
        for j in 1:size(H_r, 2)
            k_yj = k_y[j]
            for i in 1:size(H_r, 1)
                k_xi = k_x[i]
                k2 = k_xi^2 + k_yj^2
                if !(k2 ≈ zero(T))
                    phase = phase_x[i] * phase_y[j]
                    cheb_ij = cheb_mat[i, j]
                    for k in 1:size(H_r, 3)
                        r_zk = r_z[k]
                        H_r[i, j, k] += q * phase * cheb_ij(abs(z - r_zk))
                    end
                end
            end
        end
    end

    return H_r
end


@inbounds function interpolate_thin_single!(q::T, pos::NTuple{3, T}, pad_grid::Array{Complex{T}, 3}, gridinfo::GridInfo{2, T}, cheb_value::Vector{Array{T, 1}}, cheb_coefs::NTuple{2, ChebCoef{T}}, r_z::Vector{T}, Taylor_order::Int, L_z::T) where{T}

    idl = gridinfo.index_list
    pos_new = (pos[1], pos[2])

    near_id_image = image_grid_id(pos_new, gridinfo)
    near_pos_image = image_grid_pos(near_id_image, gridinfo)
    for i in 1:2
        dx = pos[i] - near_pos_image[i]
        pwcheb_eval!(dx, cheb_value[i], cheb_coefs[i])
    end

    for i in gridinfo.iter_list
        image_id = near_id_image.id .+ i
        for k in 1:size(pad_grid, 3)
            qn = q * ((pos[3] - r_z[k]) / L_z)^(2 * (Taylor_order - 1))
            pad_grid[idl[1][image_id[1]], idl[2][image_id[2]], k] += Complex{T}(qn * prod(cheb_value[j][i[j] + gridinfo.w[j] + 1] for j in 1:2))
        end
    end

    return nothing
end

@inbounds function interpolate_thin!(qs::Vector{T}, poses::Vector{NTuple{3, T}}, pad_grids::Vector{Array{Complex{T}, 3}}, gridinfo::GridInfo{2, T}, cheb_value::Vector{Array{T, 1}}, cheb_coefs::NTuple{2, ChebCoef{T}}, r_z::Vector{T}, L_z::T) where{T}

    @assert length(qs) == length(poses)
    set_zeros!.(pad_grids)

    for Taylor_order in 1:length(pad_grids)
        for i in 1:length(qs)
            interpolate_thin_single!(qs[i], poses[i], pad_grids[Taylor_order], gridinfo, cheb_value, cheb_coefs, r_z, Taylor_order, L_z)
        end
    end

    return nothing
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
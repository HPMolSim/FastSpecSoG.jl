@inbounds function gather_nu_single(q::T, pos::NTuple{3, T}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, H_c::Array{Complex{T}, 3}) where{T}
    x, y, z = pos
    L_x, L_y, L_z = L
    ϕ = zero(Complex{T})
    
    revise_phase_pos!(phase_x, phase_y, k_x, k_y, x - T(L[1] / 2), y - T(L[2] / 2))

    for k in 1:size(H_c, 3)
        cheb_val = chebpoly(k - 1, z - L_z / T(2), L_z / T(2))
        ϕ += dot(phase_x', (@view H_c[:, :, k]), phase_y)* cheb_val
    end

    return q * ϕ / (2 * L_x * L_y)
end

function gather_nu(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, H_c::Array{Complex{T}, 3}) where{T}
    E = zero(Complex{T})
    for i in 1:length(qs)
        E += gather_nu_single(qs[i], poses[i], L, k_x, k_y, phase_x, phase_y, H_c)
    end
    return real(E)
end

@inbounds function gather_thin_single(q::T, pos::NTuple{3, T}, L_z::T, pad_grid::Array{Complex{T}, 3}, gridinfo::GridInfo{2, T}, cheb_value::Vector{Array{T, 1}}, cheb_coefs::NTuple{2, ChebCoef{T}}) where{T}
    
    potential_i = zero(T)
    idl = gridinfo.index_list
    x, y, z = pos

    near_id_image = image_grid_id((pos[1], pos[2]), gridinfo)
    near_pos_image = image_grid_pos(near_id_image, gridinfo)
    for i in 1:2
        dx = pos[i] - near_pos_image[i]
        pwcheb_eval!(dx, cheb_value[i], cheb_coefs[i])
    end

    for i in gridinfo.iter_list
        image_id = near_id_image.id .+ i
        for k in 1:size(pad_grid, 3)
            cheb_val = chebpoly(k - 1, z - L_z / T(2), L_z / T(2))
            potential_i += real(pad_grid[idl[1][image_id[1]], idl[2][image_id[2]], k]) * prod(cheb_value[j][i[j] + gridinfo.w[j] + 1] for j in 1:2) * cheb_val
        end
    end

    return q * prod(gridinfo.h) * potential_i
end

function gather_thin(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, H_c::Array{Complex{T}, 3}, gridinfo::GridInfo{2, T}, cheb_value::Vector{Array{T, 1}}, cheb_coefs::NTuple{2, ChebCoef{T}}) where{T}
    E = zero(Complex{T})
    for i in 1:length(qs)
        E += gather_thin_single(qs[i], poses[i], L[3], H_c, gridinfo, cheb_value, cheb_coefs)
    end
    return real(E) / 2
end
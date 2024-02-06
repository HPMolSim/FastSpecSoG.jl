@inbounds function interpolate_nu_single!(q::T, pos::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, r_z::Vector{T}, us_mat::Array{T, 3}, H_r::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    x, y, z = pos
    for i in 1:size(H_r, 1), j in 1:size(H_r, 2), k in 1:size(H_r, 3)
        k_xi = k_x[i]
        k_yj = k_y[j]
        r_zk = r_z[k]

        phase = exp( - T(1)im * (k_xi * x + k_yj * y))

        val = zero(T)

        for l in M_mid + 1:length(uspara.sw)
            sl, wl = uspara.sw[l]
            val += wl * (T(2) - T(4) * (z - r_zk)^2 / sl^2 + (k_xi^2 + k_yj^2) *  sl^2) * exp(- (z - r_zk)^2 / sl^2) * us_mat[i, j, l - M_mid]
        end

        H_r[i, j, k] += q * π * phase * val
    end 

    return H_r
end

function interpolate_nu!(qs::Vector{T}, poses::Vector{NTuple{3, T}}, k_x::Vector{T}, k_y::Vector{T}, r_z::Vector{T}, us_mat::Array{T, 3}, H_r::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    set_zeros!(H_r)
    for i in 1:length(qs)
        interpolate_nu_single!(qs[i], poses[i], k_x, k_y, r_z, us_mat, H_r, uspara, M_mid)
    end

    return H_r
end

@inbounds function real2Cheb!(H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

    set_zeros!(H_c)
    N_z = size(H_r, 3)
    for i in 1:size(H_r, 1), j in 1:size(H_r, 2), k in 1:N_z, l in 1:N_z
        H_c[i, j, k] += 2 / N_z * H_r[i, j, l] * chebpoly(k - 1, r_z[l] - L_z / T(2), L_z / T(2))
        # cos((k - 1)* acos(r_z[l] / (L_z / 2) - 1.0))
    end

    return H_c
end
@inline function Fourier_k(i::Int, N::Int, L::T) where{T}
    return T(2π * (i - 1 - N) / L)
end

@inbounds function interpolate_nu_single!(q::T, pos::NTuple{3, T}, N::NTuple{3, Int}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, r_z::Vector{T}, us_mat::Array{T, 3}, H_r::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    x, y, z = pos
    revise_phase_neg!(phase_x, phase_y, k_x, k_y, x, y)

    for k in 1:size(H_r, 3)
        r_zk = r_z[k]
        for l in M_mid + 1:length(uspara.sw)
            sl, wl = uspara.sw[l]
            exp_temp = exp(- (z - r_zk)^2 / sl^2)
            z_temp = T(2) - T(4) * (z - r_zk)^2 / sl^2

            for j in 1:size(H_r, 2)
                k_yj = Fourier_k(j, N[2], L[2])
                for i in 1:size(H_r, 1)
                    k_xi = Fourier_k(i, N[1], L[1])
                    phase = phase_x[i] * phase_y[j]
                    us = us_mat[i, j, l - M_mid]
                    H_r[i, j, k] += q * π * phase * (z_temp + (k_xi^2 + k_yj^2) * sl^2) * exp_temp * us
                end
            end
        end
    end

    return H_r
end

function interpolate_nu!(qs::Vector{T}, poses::Vector{NTuple{3, T}}, N::NTuple{3, Int}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, r_z::Vector{T}, us_mat::Array{T, 3}, H_r::Array{Complex{T}, 3}, uspara::USeriesPara{T}, M_mid::Int) where{T}

    set_zeros!(H_r)
    for i in 1:length(qs)
        interpolate_nu_single!(qs[i], poses[i], N, L, k_x, k_y, phase_x, phase_y, r_z, us_mat, H_r, uspara, M_mid)
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
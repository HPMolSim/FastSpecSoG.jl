@inbounds function boundaries_single!(q::T, pos::NTuple{3, T}, b_l::Array{Complex{T}, 2}, b_u::Array{Complex{T}, 2}, k_x::Vector{T}, k_y::Vector{T}, phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, L_z::T, uspara::USeriesPara{T}, M_mid::Int) where{T}

    x, y, z = pos

    revise_phase_neg!(phase_x, phase_y, k_x, k_y, x, y)

    for l in M_mid + 1:length(uspara.sw)
        sl, wl = uspara.sw[l]
        temp = q * π * wl * sl^2
        temp_l = temp * exp(- z^2  / sl^2)
        temp_u = temp * exp(- (L_z - z)^2  / sl^2)
        for j in size(b_l, 2)
            k_yj = k_y[j]
            for i in 1:size(b_l, 1)
                k_xi = k_x[i]
                phase = phase_x[i] * phase_y[j]
                exp_temp = exp(-sl^2 * (k_xi^2 + k_yj^2) / 4)
                b_l[i, j] += temp_l * exp_temp * phase
                b_u[i, j] += temp_u * exp_temp * phase
            end
        end
    end

    return b_l, b_u
end

function boundaries!(qs::Vector{T}, poses::Vector{NTuple{3, T}}, b_l::Array{Complex{T}, 2}, b_u::Array{Complex{T}, 2}, k_x::Vector{T}, k_y::Vector{T}, phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, L_z::T, uspara::USeriesPara{T}, M_mid::Int) where{T}

    set_zeros!(b_l)
    set_zeros!(b_u)

    for i in 1:length(qs)
        boundaries_single!(qs[i], poses[i], b_l, b_u, k_x, k_y, phase_x, phase_y, L_z, uspara, M_mid)
    end

    return b_l, b_u
end

@inbounds function solve_eqs!(rhs::Vector{Complex{T}}, sol::Vector{Complex{T}}, H_c::Array{Complex{T}, 3}, H_s::Array{Complex{T}, 3}, b_l::Array{Complex{T}, 2}, b_u::Array{Complex{T}, 2}, ivsm::Array{T, 4}, L_z::T) where{T}

    set_zeros!(H_s)

    N_x = Int((size(H_c, 1) - 1)/2)
    N_y = Int((size(H_c, 2) - 1)/2)
    N_z = size(H_c, 3)
    for j in 1:size(H_c, 2), i in 1:size(H_c, 1)
        if !((i == N_x + 1) && (j == N_y + 1))
            for l in 1:N_z
                rhs[l] = - H_c[i, j, l]
            end
            rhs[N_z + 1] = b_l[i, j]
            rhs[N_z + 2] = b_u[i, j]

            mul!(sol, (@view ivsm[:, :, i, j]), rhs, true, false)

            H_s[i, j, 1] = T(2) * sol[N_z + 1]
            H_s[i, j, 2] = L_z^2 * (sol[2] - sol[4]) / T(32) + L_z * sol[N_z + 2] / T(2)

            sol[N_z + 1] = zero(Complex{T})
            sol[N_z + 2] = zero(Complex{T})

            for s in 3:N_z
                k = s - 1
                H_s[i, j, s] = L_z^2 / (8 * k) * 
                ((sol[s - 2] - sol[s])/(2 * (k - 1)) - (sol[s] - sol[s + 2])/(2 * (k + 1)))
            end
        end
    end

    return H_s
end
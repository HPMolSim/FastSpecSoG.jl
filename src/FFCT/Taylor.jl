function Greens_Q2d_direct(z::T, k_xi::T, k_yj::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    t = zero(T)
    k2 = k_xi^2 + k_yj^2

    for l in M_mid + 1:length(uspara.sw)
        sl, wl = uspara.sw[l]
        t += π * wl * sl^2 * exp(-sl^2 * k2 / 4) * exp(-z^2 / sl^2)
    end

    return t
end


#(x^n / n!)
@inbounds @inline function special_powern(x::T, n::Int) where{T}
    t = one(T)
    for i in 1:n
        t *= x / T(i)
    end
    return t
end

# the resulting Polynomial is in the form of f = x -> poly((x / L_z)^2)
function TaylorUSeries_k(k_xi::T, k_yj::T, L_z::T, uspara::USeriesPara{T}, M_mid::Int, Q::Int) where{T}
    coefs = zeros(T, Q)
    k2 = k_xi^2 + k_yj^2
    if !(k2 ≈ zero(T))
        for l in M_mid + 1:length(uspara.sw)
            sl, wl = uspara.sw[l]
            for n in 0:Q - 1
                coefs[n + 1] += (-1)^n * π * wl * sl^2 * exp(-sl^2 * k2 / 4) * special_powern(L_z^2 / sl^2, n)
            end
        end
    end
    return coefs
end

function TaylorUSeries(k_x::Vector{T}, k_y::Vector{T}, L_z::T, uspara::USeriesPara{T}, M_mid::Int, Q::Int) where{T}

    taylor_mats = [zeros(T, (length(k_x), length(k_y))) for i in 1:Q]

    for i in 1:length(k_x)
        for j in 1:length(k_y)
            coefs = TaylorUSeries_k(k_x[i], k_y[j], L_z, uspara, M_mid, Q)
            for k in 1:Q
                taylor_mats[k][i, j] = coefs[k]
            end
        end
    end

    return taylor_mats
end
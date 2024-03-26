function USeries_direct(z::T, k_xi::T, k_yj::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    t = zero(T)
    k2 = k_xi^2 + k_yj^2

    for l in M_mid + 1:length(uspara.sw)
        sl, wl = uspara.sw[l]
        t += π * wl * exp(-sl^2 * k2 / 4) * exp(-z^2 / sl^2) * (T(2) - T(4) * z^2 / sl^2 + k2 * sl^2)
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

function TaylorUSeries_k(k_xi::T, k_yj::T, L_z::T, uspara::USeriesPara{T}, M_mid::Int, Q::Int) where{T}
    coefs = zeros(T, Q + 1)
    k2 = k_xi^2 + k_yj^2
    for l in M_mid + 1:length(uspara.sw)
        sl, wl = uspara.sw[l]
        for n in 0:Q - 1
            coefs[n + 1] += (-1)^n * wl / 2 * (T(2) + 4n + k2 * sl^2) * exp(-sl^2 * k2 / 4) * special_powern(L_z^2 / sl^2, n)
        end
        coefs[Q + 1] += (-1)^Q * wl * Q * exp(-sl^2 * k2 / 4) * special_powern(L_z^2 / sl^2, Q)
    end
    return Polynomial(coefs)
end

function TaylorUSeries(k_x::Vector{T}, k_y::Vector{T}, L_z::T, uspara::USeriesPara{T}, M_mid::Int, Q::Int) where{T}

    taylor_mat = Array{Polynomial{T}, 2}(undef, (length(k_x), length(k_y)))

    for i in 1:length(k_x)
        for j in 1:length(k_y)
            taylor_mat[i, j] = TaylorUSeries_k(k_x[i], k_y[j], L_z, uspara, M_mid, Q)
        end
    end

    return taylor_mat
end
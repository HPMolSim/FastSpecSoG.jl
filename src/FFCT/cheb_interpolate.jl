function USeries_direct(z::T, k_xi::T, k_yj::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    t = zero(T)
    k2 = k_xi^2 + k_yj^2

    for l in M_mid + 1:length(uspara.sw)
        sl, wl = uspara.sw[l]
        t += π * wl * exp(-sl^2 * k2 / 4) * exp(-z^2 / sl^2) * (T(2) - T(4) * z^2 / sl^2 + k2 * sl^2)
    end

    return t
end

function ChebUSeries_k(k_xi::T, k_yj::T, L_z::T, uspara::USeriesPara{T}, M::Int, Q::Int) where{T}

    f = z -> USeries_direct(z, k_xi, k_yj, uspara, M)
    x = chebpoints(Q, zero(T), L_z)

    return chebinterp(f.(x), zero(T), L_z)
end

function ChebUSeries(k_x::Vector{T}, k_y::Vector{T}, L_z::T, uspara::USeriesPara{T}, M::Int, Q::Int) where{T}

    cheb_mat = Array{ChebPoly{1, T, T}, 2}(undef, (length(k_x), length(k_y)))

    for i in 1:length(k_x)
        for j in 1:length(k_y)
            cheb_mat[i, j] = ChebUSeries_k(k_x[i], k_y[j], L_z, uspara, M, Q)
        end
    end

    return cheb_mat
end
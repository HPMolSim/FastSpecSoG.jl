function USeries_direct_0(z::T, L_z::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    t = zero(T)

    for l in M_mid + 1:length(uspara.sw)
        sl, wl = uspara.sw[l]
        N = proper_N(1e-16, (L_z / sl)^2)
        
        if N == 0
            t += π * wl * sl^2 * exp(-z^2 / sl^2)
        else
            for k in 1:N
                t += π * wl * sl^2 * special_powern(- z^2 / sl^2, k)
            end
        end
    end

    return t
end

function USeries_direct(z::T, k_xi::T, k_yj::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    t = zero(T)
    k2 = k_xi^2 + k_yj^2

    for l in M_mid + 1:length(uspara.sw)
        sl, wl = uspara.sw[l]
        t += π * wl * sl^2 * exp(-sl^2 * k2 / 4) * exp(-z^2 / sl^2)
    end

    return t
end

function ChebUSeries_0(L_z::T, uspara::USeriesPara{T}, M::Int, Q::Int) where{T}

    f = z -> USeries_direct_0(z, L_z, uspara, M)
    x = chebpoints(Q, zero(T), L_z)

    return chebinterp(f.(x), zero(T), L_z)
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

@inline function chebpoly(n::Int, x::T, scale::T) where{T}
    return cos(n * acos(x / scale))
end

@inbounds function real2Cheb!(H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

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

@inbounds function Cheb2real!(H_c::Array{Complex{T}, 3}, H_r::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

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
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

@inline function chebpoly(n::Int, x::T, scale::T) where{T}
    return cos(n * acos(x / scale))
end

function dct(x::Vector{T}, fx::Vector{T}, scale::T) where{T}
    N = length(x)
    a = zeros(T, N)
    for m in 0:N-1
        for n in 1:N
            cheb_temp = m == 0 ? 0.5 : cos(m * acos(x[n] / scale))
            a[m + 1] += 2 / N * fx[n] * cheb_temp
        end
    end
    return a
end

function idct(a::Vector{T}, x::T, scale::T) where{T}
    N = length(a)
    y = zero(T)
    for m in 0:N-1
        y += a[m + 1] * cos(m * acos(x / scale))
    end
    return y
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

@inbounds function real2Cheb_Q2D!(H_r::Array{Complex{T}, 3}, H_c::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

    set_zeros!(H_c)
    N_z = size(H_r, 3)

    # for k in 1:N_z, l in 1:N_z
    #     cheb_temp = k == 1 ? 0.5 : chebpoly(k - 1, r_z[l] - L_z / T(2), L_z / T(2))
    #     for j in 1:size(H_r, 2), i in 1:size(H_r, 1)
    #         H_c[i, j, k] += 2 / N_z * H_r[i, j, l] * cheb_temp
    #     end
    # end

    for j in 1:size(H_r, 2), i in 1:size(H_r, 1)
        (@view H_c[i, j, :]) .= chebinterp((@view H_c[i, j, :]), zero(T), L_z).coefs
    end

    return H_c
end

@inbounds function Cheb2real_Q2D!(H_c::Array{Complex{T}, 3}, H_r::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

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
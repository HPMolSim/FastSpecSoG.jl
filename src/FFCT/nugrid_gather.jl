@inbounds function gather_nu_single(q::T, pos::NTuple{3, T}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, H_s::Array{Complex{T}, 3}) where{T}
    x, y, z = pos
    L_x, L_y, L_z = L
    ϕ = zero(Complex{T})
    
    
    for i in 1:size(H_s, 1), j in 1:size(H_s, 2)
        k_xi = k_x[i]
        k_yj = k_y[j]
        phase = exp(T(1)im * (k_xi * x + k_yj * y))

        ϕ += H_s[i, j, 1] * phase / T(2)
        for k in 2:size(H_s, 3)
            ϕ += H_s[i, j, k] * phase * chebpoly(k - 1, z - L_z / T(2), L_z / T(2))
        end
    end

    return q * ϕ / (2 * L_x * L_y)
end

function gather_nu(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, H_s::Array{Complex{T}, 3}) where{T}
    E = zero(Complex{T})
    for i in 1:length(qs)
        E += gather_nu_single(qs[i], poses[i], L, k_x, k_y, H_s)
    end
    return real(E)
end
function BSAj(x::T, b::T, σ::T, j::Int) where{T}
    return T(2.0) * log(b) / sqrt(T(2π)) / σ / b^j * exp(- (x / b^j / σ)^2 / T(2.0))
end

function BSA(x::T, b::T, σ::T, M::Int) where{T}
    BSA_value = zero(T)
    for j in -M:M
        BSA_value += BSAj(x, b, σ, j)
    end
    return BSA_value
end

function U_series(x::T, b::T, σ::T, M::Int) where{T}
    U_series_value = zero(T)

    for l in 0:M
        U_series_value += exp(- x^2 / (2 * b^(2 * l) * σ^2)) / b^l
    end

    U_series_value *= sqrt(2 / π) * log(b) / σ

    return U_series_value
end

function USeriesPara(b::T, σ::T, ω::T, M::Int) where{T}
    sw = Vector{NTuple{2, T}}(undef, M + 1)

    sw[1] = (sqrt(2 * σ^2), ω * sqrt(2 / π) * log(b) / σ)

    for l in 1:M
        sw[l + 1] = (sqrt(2 * b^(2 * l) * σ^2), sqrt(2 / π) * log(b) / σ / b^l)
    end
    return USeriesPara(sw)
end

function U_series(x::T, uspara::USeriesPara{T}) where{T}
    U_series_value = zero(T)

    for i in eachindex(uspara.sw)
        (s, w) = uspara.sw[i]
        U_series_value += w * exp(- x^2 / s^2)
    end

    return U_series_value
end
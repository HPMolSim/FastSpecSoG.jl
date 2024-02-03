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

const preset_parameters = [
    tuple(2.0, 5.027010924194599, 0.994446492762232252, 6), 
    tuple(1.62976708826776469, 3.633717409009413, 1.00780697934380681, 16), 
    tuple(1.48783512395703226, 2.662784519725113, 0.991911705759818, 30), 
    tuple(1.32070036405934420, 2.277149356440992, 1.00188914114811980, 64),
    tuple(1.21812525709410644, 1.774456369233284, 1.00090146156033341, 102)]

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

function USeriesPara(preset::Int, T::DataType = Float64)
    @assert preset ≤ length(preset_parameters)

    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])

    return USeriesPara(b, σ, ω, M)
end

function U_series(x::T, uspara::USeriesPara{T}) where{T}
    U_series_value = zero(T)

    for i in eachindex(uspara.sw)
        (s, w) = uspara.sw[i]
        U_series_value += w * exp(- x^2 / s^2)
    end

    return U_series_value
end
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
    tuple(2.0, 0.5027010924194599, 0.994446492762232252, 18), 
    tuple(1.62976708826776469, 0.3633717409009413, 1.00780697934380681, 39), 
    tuple(1.48783512395703226, 0.2662784519725113, 0.991911705759818, 45), 
    tuple(1.32070036405934420, 0.2277149356440992, 1.00188914114811980, 78),
    tuple(1.21812525709410644, 0.1774456369233284, 1.00090146156033341, 144),
    tuple(1.14878150173321925, 0.1370684570284264 , 1.000036834835822, 233)]

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

function USeriesPara(preset::Int; r_c::T = 10.0) where{T}
    @assert preset ≤ length(preset_parameters)

    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]) * r_c, T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])

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

function proper_M(η::T, L_z::T, uspara::USeriesPara{T}) where{T}
    if uspara.sw[1][1] > η * L_z
        return 0
    end
    for i in 1:length(uspara.sw) - 1
        (s, w) = uspara.sw[i]
        (sp, wp) = uspara.sw[i + 1]
        if  s ≤ η * L_z ≤ sp
            return i
        end
    end
    @error "No proper M found for η = $η"
end
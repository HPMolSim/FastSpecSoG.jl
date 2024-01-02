struct USeriesPara{T}
    sw::Vector{NTuple{2, T}}
end

struct FastSpecSOGInteraction{T} <: ExTinyMD.AbstractInteraction
    b::T
    σ::T
    ω::T
    M::Int
    r_c::T
    k_c::T
    ϵ::T
    L::NTuple{3, T}
    k_set::Vector{NTuple{3, T}}
    uspara::USeriesPara{T}
    n_atoms::Int
end

function k_set_2D(k_c::T, L::NTuple{3,T}) where{T}
    L_x, L_y, L_z = L

    mx_max = ceil(Int, k_c * L_x / 2π) + 1
    my_max = ceil(Int, k_c * L_y / 2π) + 1

    k_set = Vector{NTuple{3, T}}()

    for m_x in - mx_max : mx_max
        for m_y in - my_max : my_max
            k_x = m_x * 2π / L_x
            k_y = m_y * 2π / L_y
            k = sqrt(k_x^2 + k_y^2)
            if k <= k_c
                push!(k_set, (k_x, k_y, k))
            end
        end
    end

    return k_set
end


function FastSpecSOGInteraction(L::NTuple{3, T}, n_atoms::Int, r_c::T, k_c::T, b::T, σ::T, ω::T, M::T; ϵ::T = 1.0) where{T}
    uspara = USeriesPara(b, σ, ω, M)
    k_set = k_set_2D(k_c, L)
    return FastSpecSOGInteraction{T}(b, σ, ω, M, r_c, k_c, ϵ, L, k_set, uspara, n_atoms)
end

function FastSpecSOGInteraction(L::NTuple{3, T}, n_atoms::Int, r_c::T, k_c::T; preset::Int = 1, ϵ::T = 1.0) where{T}

    # these are preset parameters
    preset_parameters = [
        [2.0, 5.027010924194599, 0.994446492762232252, 6], 
        [1.62976708826776469, 3.633717409009413, 1.00780697934380681, 16], 
        [1.48783512395703226, 2.662784519725113, 0.991911705759818, 30], 
        [1.32070036405934420, 2.277149356440992, 1.00188914114811980, 64],
        [1.21812525709410644, 1.774456369233284, 1.00090146156033341, 102]]

    @assert preset ≤ length(preset_parameters)
    
    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])
    uspara = USeriesPara(b, σ, ω, M)
    k_set = k_set_2D(k_c, L)
    return FastSpecSOGInteraction{T}(b, σ, ω, M, r_c, k_c, ϵ, L, k_set, uspara, n_atoms)
end
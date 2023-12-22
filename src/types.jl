struct USeriesPara{T}
    sw::Vector{NTuple{2, T}}
end

struct FastSpecSOGInteraction{T} <: ExTinyMD.AbstractInteraction
    b::T
    σ::T
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


function FastSpecSOGInteraction(L::NTuple{3, T}, n_atoms::Int;ϵ::T = 1.0, b::T = 1.62976708826776469, σ::T = 3.633717409009413, M::Int = 25, r_c::T = 10.0, k_c::T = 2.0) where{T}
    uspara = USeriesPara(b, σ, M)
    k_set = k_set_2D(k_c, L)
    return FastSpecSOGInteraction{T}(b, σ, M, r_c, k_c, ϵ, L, k_set, uspara, n_atoms)
end
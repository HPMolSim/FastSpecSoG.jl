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


function FSSoG_naive(L::NTuple{3, T}, n_atoms::Int, r_c::T, k_c::T, b::T, σ::T, ω::T, M::T; ϵ::T = 1.0) where{T}
    uspara = USeriesPara(b, σ, ω, M)
    k_set = k_set_2D(k_c, L)
    return FSSoG_naive{T}(b, σ, ω, M, r_c, k_c, ϵ, L, k_set, uspara, n_atoms)
end

function FSSoG_naive(L::NTuple{3, T}, n_atoms::Int, r_c::T, k_c::T; preset::Int = 1, ϵ::T = 1.0) where{T}

    @assert preset ≤ length(preset_parameters)
    
    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])
    uspara = USeriesPara(b, σ, ω, M)
    k_set = k_set_2D(k_c, L)
    return FSSoG_naive{T}(b, σ, ω, M, r_c, k_c, ϵ, L, k_set, uspara, n_atoms)
end

function FSSoGInteraction(
    L::NTuple{3, T}, n_atoms::Int, 
    n_k::NTuple{3, Int}, h::T, win_cut::Int, win_width::T, β::T,; 
    preset::Int = 1, ϵ::T = one(T)) where{T}
    @assert preset ≤ length(preset_parameters)
    
    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])
    uspara = USeriesPara(b, σ, ω, M)

    

    return 
end
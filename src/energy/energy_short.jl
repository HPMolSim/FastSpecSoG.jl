function F0_cal(b::T, σ::T, ω::T, M::Int) where{T}
    return - log(b) / sqrt(2π) / σ * (ω + (one(T) - b^(-M)) / (b - one(T)))
end

function F0_cal(;preset::Int = 1, T::DataType = Float64)
    @assert preset ≤ length(preset_parameters)

    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]), T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])

    return F0_cal(b, σ, ω, M)
end

function Es_USeries_Cheb(uspara::USeriesPara{T}, r_min::T, r_max::T, Q::Int) where{T}
    f = r -> U_series(r, uspara)
    x = chebpoints(Q, r_min, r_max)

    return chebinterp(f.(x), r_min, r_max)
end

function Es_Cheb_precompute(preset::Int, r_min::T, r_max::T, Q::Int) where{T}
    uspara = USeriesPara(preset, T)
    uspara_cheb = Es_USeries_Cheb(uspara, r_min, r_max, Q)
    F0 = F0_cal(preset = preset, T = T)
    return uspara_cheb, F0
end

function Es_self(q::Vector{T}, F0::T) where{T}

    Q = sum(qi^2 for qi in q)

    return Q * F0
end

function Es_Cheb_pair(q_1::T, q_2::T, uspara_cheb::ChebPoly{1, T, T}, r::T) where{T}
    return q_1 * q_2 * (one(T) / r - uspara_cheb(r))
end

function short_energy_Cheb(uspara_cheb::ChebPoly{1, T, T}, r_c::T, F0::T, boundary::Boundary{T}, neighbor::CellList3D{T}, position::Vector{NTuple{3, T}}, q::Vector{T}) where{T}
    neighbor_list = neighbor.neighbor_list

    energy_short = zero(T)

    for (i, j, ρ) in neighbor_list
        coord_1, coord_2, r_sq = position_check3D(position[i], position[j], boundary, r_c)
        if iszero(r_sq)
            nothing
        else
            q_1 = q[i]
            q_2 = q[j]
            energy_short += Es_Cheb_pair(q_1, q_2, uspara_cheb, sqrt(r_sq))
        end
    end

    energy_short += Es_self(q, F0)

    return energy_short / 4π
end

function energy_short(interaction::FSSoGInteraction{T}, neighbor::CellList3D{T}) where{T}
    return short_energy_Cheb(interaction.uspara_cheb, interaction.r_c, interaction.F0, interaction.boundary, neighbor, interaction.position, interaction.charge) / interaction.ϵ
end

function energy_short(interaction::FSSoGThinInteraction{T}, neighbor::CellList3D{T}) where{T}
    return short_energy_Cheb(interaction.uspara_cheb, interaction.r_c, interaction.F0, interaction.boundary, neighbor, interaction.position, interaction.charge) / interaction.ϵ
end
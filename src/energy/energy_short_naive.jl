function short_energy_naive(interaction::FSSoG_naive{T}, neighbor::CellList3D{T}, position::Vector{Point{3, T}}, q::Vector{T}) where{T}
    neighbor_list = neighbor.neighbor_list

    energy_short = zero(T)
    boundary = Q2dBoundary(interaction.L[1], interaction.L[2], interaction.L[3])

    for (i, j, ρ) in neighbor_list
        coord_1, coord_2, r_sq = position_check3D(position[i], position[j], boundary, interaction.r_c)
        if iszero(r_sq)
            nothing
        else
            q_1 = q[i]
            q_2 = q[j]
            energy_short += Es_naive_pair(q_1, q_2, interaction.uspara, r_sq)
        end
    end

    energy_short += Es_naive_self(q, interaction)

    return energy_short / (4π * interaction.ϵ)
end

function Es_naive_pair(q_1::T, q_2::T, uspara::USeriesPara{T}, r_sq::T) where{T}
    return q_1 * q_2 * (one(T) / sqrt(r_sq) - U_series(sqrt(r_sq), uspara))
end

function Es_naive_self(q::Vector{T}, interaction::FSSoG_naive{T}) where{T}
    Q = sum(qi^2 for qi in q)

    F0 = - log(interaction.b) / sqrt(2π) / interaction.σ * (interaction.ω + (one(T) - interaction.b^(-interaction.M)) / (interaction.b - one(T)))

    return Q * F0
end
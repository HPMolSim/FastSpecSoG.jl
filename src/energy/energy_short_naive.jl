# Kept as its own function so that each candidate-iterator type specializes.
function _short_naive_sum(pairs, interaction::FSSoG_naive{T}, position, q::Vector{T}) where{T}
    r_c_sq = interaction.r_c^2
    energy_short = zero(T)

    for pair in pairs
        i, j = pair[1], pair[2]
        _, _, r_sq = _min_image_slab(position[i], position[j], interaction.L)
        # See `_short_Cheb_sum`: the cutoff test and the coincident-pair test
        # were previously conflated in `position_check3D`'s all-zero sentinel.
        # `Es_naive_pair` computes 1/sqrt(r_sq), so r_sq == 0 must be skipped.
        (r_sq < r_c_sq && !iszero(r_sq)) || continue
        energy_short += Es_naive_pair(q[i], q[j], interaction.uspara, r_sq)
    end

    return energy_short
end

"""
    short_energy_naive(interaction, position, q; neighbor_list = nothing)

Direct (unapproximated) short-range energy, the reference the Chebyshev form in
[`short_energy_Cheb`](@ref) is validated against. `position` is
array-of-structs and only indexed.

**`interaction.r_c` must be strictly less than `min(Lx, Ly) / 2`.**
`neighbor_list` is an optional iterable of candidate pairs; omitting it falls
back to all i < j pairs.
"""
function short_energy_naive(interaction::FSSoG_naive{T}, position, q::Vector{T}; neighbor_list = nothing) where{T}

    energy_short = _short_naive_sum(_candidate_pairs(_pair_list(neighbor_list), length(q)), interaction, position, q)

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

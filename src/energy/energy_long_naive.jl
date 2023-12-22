function long_energy_naive(interaction::FastSpecSOGInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T}
    energy = zero(T)
    for K in interaction.k_set
        energy += long_energy_naive_k(K, interaction, position, charge)
    end
    return energy / (4π * interaction.ϵ)
end

function long_energy_naive_k(K::NTuple{3, T}, interaction::FastSpecSOGInteraction{T}, position::Vector{Point{3, T}}, charge::Vector{T}) where{T}

    energy_k = zero(T)

    kx, ky, k = K
    for i in 1:interaction.n_atoms
        qi = charge[i]
        xi, yi, zi = position[i]

        for j in 1:interaction.n_atoms
            ϕ_ij = zero(T)
            qj = charge[j]
            xj, yj, zj = position[j]
            for (s, w) in interaction.uspara.sw
                ϕ_ij += w * s^2 * exp( - (zi - zj)^2 / s^2) * exp(-s^2 * k^2 / 4) * cos(kx * (xi - xj) + ky * (yi - yj))
            end
            energy_k += qi * qj * ϕ_ij
        end
    end

    return energy_k * π / (interaction.L[1] * interaction.L[2])
end
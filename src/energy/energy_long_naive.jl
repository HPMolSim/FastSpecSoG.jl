"""
    long_energy_naive(interaction, position, charge)

Compute the long-range energy using a naive approach.

# Arguments
- `interaction`: The FSSoG_naive object representing the interaction parameters.
- `position`: A vector of 3D positions of the charges.
- `charge`: A vector of charges.

# Returns
The computed long-range energy.

"""
function long_energy_naive(interaction::FSSoG_naive{T}, position::Vector{NTuple{3, T}}, charge::Vector{T}) where{T}
    energy = zero(T)
    for K in interaction.k_set
        energy += long_energy_naive_k(K, interaction, position, charge)
    end
    return energy / (4π * interaction.ϵ)
end

function long_energy_naive_k(K::NTuple{3, T}, interaction::FSSoG_naive{T}, position::Vector{NTuple{3, T}}, charge::Vector{T}) where{T}

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

    return energy_k * π / (2 * interaction.L[1] * interaction.L[2])
end

function long_energy_sw_k(qs::Vector{T}, poses::Vector{NTuple{3, T}}, cutoff::Int, L::NTuple{3, T}, s::T, w::T) where{T}

    n_atoms = length(qs)

    E = Atomic{T}(zero(T))
    @threads for n in CartesianIndices((n_atoms, n_atoms))
        i, j = Tuple(n)
        qi = qs[i]
        xi, yi, zi = poses[i]
        qj = qs[j]
        xj, yj, zj = poses[j]

        ϕ_ij = zero(T)
        for m_x in -cutoff:cutoff
            kx = 2π * m_x / L[1]
            for m_y in -cutoff:cutoff
                ky = 2π * m_y / L[2]
                k = sqrt(kx^2 + ky^2)
                if !((m_x == 0) && (m_y == 0))
                    ϕ_ij += w * s^2 * exp( - (zi - zj)^2 / s^2) * exp(-s^2 * k^2 / 4) * cos(kx * (xi - xj) + ky * (yi - yj))
                end
            end
        end

        atomic_add!(E, qi * qj * ϕ_ij)
    end

    return E[] * π / (2 * L[1] * L[2])
end

function long_energy_sw_0(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, s::T, w::T) where{T}

    n_atoms = length(qs)

    E = Atomic{T}(zero(T))

    @threads for i in 1:n_atoms
        qi = qs[i]
        xi, yi, zi = poses[i]
        t = Double64(zero(T))
        for j in 1:n_atoms
            qj = qs[j]
            xj, yj, zj = poses[j]

            t += qj * exp( - Double64((zi - zj)^2 / s^2))
        end

        atomic_add!(E, qi * T(t))
    end

    return w * s^2 * E[] * π / (2 * L[1] * L[2])
end

function long_energy_us_k(qs::Vector{T}, poses::Vector{NTuple{3, T}}, cutoff::Int, L::NTuple{3, T}, uspara::USeriesPara{T}, M_min::Int, M_max::Int) where{T}
    @assert M_min ≥ 1
    @assert M_max ≤ length(uspara.sw)

    Ek = zero(T)
    for l in M_min:M_max
        s, w = uspara.sw[l]
        Ek += long_energy_sw_k(qs, poses, cutoff, L, s, w)
    end
    @debug "long range energy, direct su, k" Ek

    return Ek
end

function long_energy_us_k(qs::Vector{T}, poses::Vector{NTuple{3, T}}, accuracy::T, L::NTuple{3, T}, uspara::USeriesPara{T}, M_min::Int, M_max::Int) where{T}
    @assert M_min ≥ 1
    @assert M_max ≤ length(uspara.sw)

    Ek = zero(T)
    for l in M_min:M_max
        s, w = uspara.sw[l]
        # accuracy = exp(-k^2 * s^2 / 4)
        km = sqrt(-4 * log(accuracy) / s^2)
        # km = π * n / max{L_x, L_y, L_z}
        cutoff = ceil(Int, km * maximum(L) / 2π) + 1
        Ek += long_energy_sw_k(qs, poses, cutoff, L, s, w)
    end
    @debug "long range energy, direct su, k" Ek

    return Ek
end

function long_energy_us_0(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, uspara::USeriesPara{T}, M_min::Int, M_max::Int) where{T}
    @assert M_min ≥ 1
    @assert M_max ≤ length(uspara.sw)

    E0 = zero(T)
    for l in M_min:M_max
        s, w = uspara.sw[l]
        E0 += long_energy_sw_0(qs, poses, L, s, w)
    end
    @debug "long range energy, direct sum, zeroth mode" E0

    return E0
end
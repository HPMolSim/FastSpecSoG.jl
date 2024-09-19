export long_energy_naive_Q1D_k_sw, long_energy_naive_Q1D_k

function long_energy_naive_Q1D_k_sw(k::T, qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, s::T, w::T) where{T}

    energy_k = zero(T)
    n_atoms = length(qs)
    for i in 1:n_atoms
        qi = qs[i]
        xi, yi, zi = poses[i]

        for j in 1:n_atoms
            ϕ_ij = zero(T)
            qj = qs[j]
            xj, yj, zj = poses[j]
            ϕ_ij += w * s^2 * exp( - ((xi - xj)^2 + (yi - yj)^2) / s^2) * exp(-s^2 * k^2 / 4) * cos(k *(zi - zj))
            energy_k += qi * qj * ϕ_ij
        end
    end

    return energy_k * π / (2 * L[3])
end

function long_energy_naive_Q1D_k(qs::Vector{T}, poses::Vector{NTuple{3, T}}, accuracy::T, L::NTuple{3, T}, uspara::USeriesPara{T}, M_min::Int, M_max::Int) where{T}
    @assert M_min ≥ 1
    @assert M_max ≤ length(uspara.sw)

    Ek = zero(T)
    for l in M_min:M_max
        s, w = uspara.sw[l]
        km = sqrt(-4 * log(accuracy) / s^2)
        cutoff = ceil(Int, km * (L[3]) / 2π) + 1
        for m in -cutoff:cutoff
            if m != 0
                k = 2π * m / L[3]
                Ek += long_energy_naive_Q1D_k_sw(k, qs, poses, L, s, w)
            end
        end
    end

    return Ek
end
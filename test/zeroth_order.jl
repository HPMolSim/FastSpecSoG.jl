function long_energy_us_0_big(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, uspara::USeriesPara{T}, M_min::Int, M_max::Int) where{T}
    @assert M_min ≥ 1
    @assert M_max ≤ length(uspara.sw)
    
    n_atoms = length(qs)

    E0 = zero(T)
    for l in M_min:M_max
        s, w = uspara.sw[l]
        for i in 1:n_atoms
            qi = qs[i]
            xi, yi, zi = poses[i]
            t = big(zero(T))
            for j in 1:n_atoms
                qj = qs[j]
                xj, yj, zj = poses[j]
    
                t += qj * exp( - big((zi - zj)^2 / s^2))
            end
    
            E0 += w * s^2 * qi * T(t) * π / (2 * L[1] * L[2])
        end
    end

    return E0
end

@testset "zeroth_order" begin
    @info "testing the zeroth_order energy"
    n_atoms = 32
    L = 20.0
    Q_0 = 64

    r_z, grids, chebcoefs = zero_paras_gen(L, Q_0)

    qs = [(-1.0)^i for i in 1:n_atoms]
    poses = [(rand() * L, rand() * L, rand() * L) for i in 1:n_atoms]

    for preset in 1:6
        uspara = USeriesPara(preset)
        M_mid = 10
        chebuseries = ChebUSeries_0(L, uspara, M_mid, Q_0)
        E_FGT = zeroth_order(qs, poses, (L, L, L), chebuseries, r_z, chebcoefs, grids)
        E_naive = long_energy_us_0(qs, poses, (L, L, L), uspara, M_mid + 1, length(uspara.sw))
        E_naive_big = long_energy_us_0_big(qs, poses, (L, L, L), uspara, M_mid + 1, length(uspara.sw))
        # @show preset, E_FGT, E_naive, (E_FGT - E_naive) / E_naive
        @test E_FGT ≈ E_naive
        @test isapprox(E_naive, E_naive_big, atol=1e-14)
    end
end
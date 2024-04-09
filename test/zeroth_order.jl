@testset "FGT1d by mixed taylor chebshev" begin
    @info "testing the FGT1d by mixed taylor chebshev"
    n_atoms = 32
    L = 20.0
    Q_0 = 32

    r_z, grids, chebcoefs = zero_paras_gen(L, Q_0)

    qs = [(-1.0)^i for i in 1:n_atoms]
    poses = [(rand() * L, rand() * L, rand() * L) for i in 1:n_atoms]

    for preset in 1:6
        uspara = USeriesPara(preset)

        for l in 5:length(uspara.sw)
            s, w = uspara.sw[l]
            E_FGT = FGT1d(qs, poses, s, L, r_z, chebcoefs, grids)
            E_naive = FGT1d_naive(qs, poses, s)
            N = FastSpecSoG.proper_N(1e-16, (L / s)^2)
            # @show l, N, E_FGT, E_naive, (E_FGT - E_naive) / E_naive
            @test E_FGT ≈ E_naive
        end
    end
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
        E_FGT = zeroth_order(qs, poses, (L, L, L), uspara, M_mid, r_z, chebcoefs, grids)
        E_naive = long_energy_us_0(qs, poses, (L, L, L), uspara, M_mid + 1, length(uspara.sw))
        # @show preset, E_FGT, E_naive, (E_FGT - E_naive) / E_naive
        @test E_FGT ≈ E_naive
    end
end
@testset "zeroth_order by SOE" begin

    n_atoms = 100
    L = 20.0
    poses = [ (rand() * L, rand() * L, rand() * L) for i in 1:n_atoms ]
    qs = [ (-1.0)^i for i in 1:n_atoms ]
    z = [ poses[i][3] for i in 1:n_atoms ]
    sort_z = sortperm(z)
    soepara = SoePara16()

    for preset in 1:6
        uspara = USeriesPara(preset)
        M_mid = 0

        E0_soe = zeroth_order(qs, z, soepara, uspara, sort_z, (L, L, L), M_mid)
        E0_direct = long_energy_us_0(qs, poses, (L, L, L), uspara, 1, length(uspara.sw))

        @show E0_soe, E0_direct, abs(E0_soe - E0_direct)
        # @test abs(E0_soe - E0_direct) < 1e-14
    end
end
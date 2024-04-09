@testset "FGT1d by chebshev" begin
    n_atoms = 256
    L = 20.0
    Q0 = 64

    qs = [(-1.0)^i for i in 1:n_atoms]
    poses = [(rand() * L, rand() * L, rand() * L) for i in 1:n_atoms]

    uspara = USeriesPara(6)

    for l in 10:length(uspara.sw)
        s, w = uspara.sw[l]
        E_FGT = FGT1d(qs, poses, s, L, Q0)
        E_naive = FGT1d_naive(qs, poses, s)
        # @show l, E_FGT, E_naive, E_FGT - E_naive
        @test E_FGT ≈ E_naive
    end
end
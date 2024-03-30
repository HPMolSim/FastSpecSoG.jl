@testset "Taylor expansions" begin
    @info "Testing the TaylorUSeries_k"
    for preset in 1:6
        @testset "Preset $preset" begin
            uspara = USeriesPara(preset)
            M_mid = 0
            L_z = abs(1.0 + randn())
            Q = 64
            k_x, k_y = randn(2)
            poly = Polynomial(TaylorUSeries_k(k_x, k_y, L_z, uspara, M_mid, Q))
            xs = [0.0:0.01:L_z...]
            fd = x -> FastSpecSoG.Greens_Q2d_direct(x, k_x, k_y, uspara, M_mid)
            fp = x -> poly((x / L_z)^2)
            @test fp.(xs) ≈ fd.(xs)
        end
    end
end
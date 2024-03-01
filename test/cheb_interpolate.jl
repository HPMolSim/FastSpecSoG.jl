@testset "Cheb expansion of U_series" begin
    @info "testing the Cheb expansion"
    k_xi = 0.1
    k_yj = 0.2
    L_z = 10.0
    z = [0.01:0.01:L_z...]
    f_cheb = FastSpecSoG.ChebUSeries_k(k_xi, k_yj, L_z, USeriesPara(2), 3, 16)
    direct_U = map(x -> FastSpecSoG.USeries_direct(x, k_xi, k_yj, USeriesPara(2), 3), z)
    @test isapprox(f_cheb.(z), direct_U)
end
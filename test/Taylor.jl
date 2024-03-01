@testset "Taylor expansion by horner's rule" begin
    @info "testing the Taylor expansion"
    for n in 1:10
        coefs = randn(n)
        taylorseries = TaylorSeriesPara(coefs)
        f = Taylor1(coefs)
        x = randn()
        @test isapprox(horner(x, taylorseries), f(x))
    end
end

@testset "Taylor expansion of U_series" begin
    @info "testing the Taylor expansion"
    k_xi = 0.1
    k_yj = 0.2
    L_z = 10.0
    z = 1.0
    taylor_U = TaylorUSeries_k(k_xi, k_yj, L_z, USeriesPara(2), 3, 16)
    direct_U = USeries_direct(z, k_xi, k_yj, USeriesPara(2), 3)
    @test isapprox(horner((z/L_z)^2, taylor_U), direct_U)
end
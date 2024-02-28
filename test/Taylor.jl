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
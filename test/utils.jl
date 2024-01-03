function d_poly(x, p)
    r = 0.0
    for i in 0:length(p) - 1
        r += x^i * p[i + 1]
    end
    return r
end

@testset "test horner" begin
    for N in 1:10
        p = rand(N)
        x = rand()
        @test horner(x, p) ≈ d_poly(x, p)
    end
end
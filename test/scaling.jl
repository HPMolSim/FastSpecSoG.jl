@testset "fft_Q2D" begin
    N = 32
    grid = rand(ComplexF64, N, N, N)
    grid_copy = copy(grid)
    grid = piecewise_fft!(grid)
    for k in 1:N
        @test isapprox(grid[:, :, k], fft(grid_copy[:, :, k]))
    end
end

@testset "ifft_Q2D" begin
    N = 32
    grid = rand(ComplexF64, N, N, N)
    grid_copy = copy(grid)
    grid = piecewise_ifft!(grid)
    for k in 1:N
        @test isapprox(grid[:, :, k], ifft(grid_copy[:, :, k]))
    end
end

@testset "scaling_Q2D" begin
    N = 32
    grid = rand(ComplexF64, N, N, N)
    factor = rand(ComplexF64, N, N)
    grid_copy = copy(grid)
    grid = piecewise_mul!(grid, factor)
    for k in 1:N
        @test isapprox(grid[:, :, k], grid_copy[:, :, k] .* factor)
    end
end

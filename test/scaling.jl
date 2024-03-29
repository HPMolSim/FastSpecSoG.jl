@testset "fft_Q2D" begin
    N = 32
    grid = rand(ComplexF64, N, N, N)
    grid_copy = copy(grid)
    grid = fft_Q2d!(grid)
    for k in 1:N
        @test isapprox(grid[:, :, k], fft(grid_copy[:, :, k]))
    end
end

@testset "ifft_Q2D" begin
    N = 32
    grid = rand(ComplexF64, N, N, N)
    grid_copy = copy(grid)
    grid = ifft_Q2d!(grid)
    for k in 1:N
        @test isapprox(grid[:, :, k], ifft(grid_copy[:, :, k]))
    end
end

@testset "scaling_Q2D" begin
    N = 32
    grid = rand(ComplexF64, N, N, N)
    factor = rand(ComplexF64, N, N)
    grid_copy = copy(grid)
    grid = scaling_Q2D!(grid, factor)
    for k in 1:N
        @test isapprox(grid[:, :, k], grid_copy[:, :, k] .* factor)
    end
end

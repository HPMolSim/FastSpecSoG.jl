@testset "compare the BSA formula" begin
    @info "testing the BSA formula"
    for (preset, accuracy) in zip([1, 2, 3, 4, 5, 6], [1e-1, 1e-3, 1e-4, 1e-6, 1e-8, 1e-12])
        @testset "preset = $preset" begin
            b, σ, ω, M = FastSpecSoG.preset_parameters[preset]
            for x in range(1.0, 10.0, length=100)
                @test abs(BSA(x, b, σ, M) * x - 1.0) < accuracy
            end
        end
    end
end
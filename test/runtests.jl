using FastSpecSoG
using Test
using ExTinyMD, EwaldSummations

@testset "FastSpecSoG.jl" begin
    include("U_series.jl")
    include("energy_naive.jl")
    include("energy_mid.jl")
    include("energy_long.jl")
end
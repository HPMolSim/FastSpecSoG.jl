using FastSpecSoG
using Test
using ExTinyMD, EwaldSummations, ProgressBars

@testset "FastSpecSoG.jl" begin
    include("U_series.jl")
    # include("energy_naive.jl")
    include("energy_long_mid.jl")
end
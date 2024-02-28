using FastSpecSoG
using Test
using ExTinyMD, EwaldSummations, TaylorSeries

@testset "FastSpecSoG.jl" begin
    include("U_series.jl")
    # include("energy_naive.jl")
    include("energy_mid.jl")
    include("interpolate_nu.jl")
    include("energy_long.jl")
    include("Taylor.jl")
end
using FastSpecSoG
using Test
using ExTinyMD, EwaldSummations, TaylorSeries
using Random
Random.seed!(1234)

@testset "FastSpecSoG.jl" begin
    include("U_series.jl")
    include("energy.jl")
    include("energy_short.jl")
    include("energy_mid.jl")
    include("interpolate_nu.jl")
    include("energy_long.jl")
    include("cheb_interpolate.jl")
end
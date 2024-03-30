using FastSpecSoG
using Test
using ExTinyMD, EwaldSummations, TaylorSeries, FFTW, Polynomials
using Random
# Random.seed!(1234)

@testset "FastSpecSoG.jl" begin
    include("U_series.jl")
    include("energy.jl")
    include("energy_short.jl")
    include("energy_mid.jl")
    include("interpolate.jl")
    include("energy_long.jl")
    include("Chebyshev.jl")
    include("Taylor.jl")
    include("energy_thin.jl")
    include("scaling.jl")
end
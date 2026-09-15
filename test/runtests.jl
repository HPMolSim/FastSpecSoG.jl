using FastSpecSoG
using Test
using ExTinyMD, TaylorSeries, FFTW, Polynomials, SpecialFunctions
using DoubleFloats
using Random
Random.seed!(1234)

@testset "FastSpecSoG.jl" begin
    include("testutils.jl")
    include("U_series.jl")
    include("plan.jl")
    include("adapter.jl")
    include("energy.jl")
    include("energy_short.jl")
    include("energy_mid.jl")
    include("interpolate.jl")
    include("energy_long.jl")
    include("Chebyshev.jl")
    include("Taylor.jl")
    include("scaling.jl")
    include("zeroth_order.jl")
    include("standalone.jl")
end
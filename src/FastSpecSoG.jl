module FastSpecSoG

using ExTinyMD, LinearAlgebra, FFTW, LoopVectorization

export USeriesPara, U_series
export FastSpecSOGInteraction
export short_energy_naive, long_energy_naive, energy_naive

include("types.jl")
include("U_series.jl")

include("energy/energy.jl")
include("energy/energy_short_naive.jl")
include("energy/energy_long_naive.jl")

end

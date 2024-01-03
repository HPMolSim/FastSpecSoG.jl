module FastSpecSoG

using ExTinyMD, LinearAlgebra, FFTW, LoopVectorization, FastTransforms, SpecialFunctions

export USeriesPara, U_series
export FSSoG_naive
export short_energy_naive, long_energy_naive, energy_naive

include("types.jl")

include("utils/U_series.jl")
include("utils/polys.jl")
include("utils/mesh.jl")

include("utils/FSSoGInteraction.jl")


include("energy/energy.jl")
include("energy/energy_short_naive.jl")
include("energy/energy_long_naive.jl")

end

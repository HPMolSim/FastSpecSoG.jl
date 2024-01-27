module FastSpecSoG

using ExTinyMD, LinearAlgebra, SpecialFunctions, ChebParticleMesh

export USeriesPara, U_series, BSA
export FSSoG_naive
export short_energy_naive, long_energy_naive, energy_naive

include("types.jl")

include("U_series.jl")
include("FSSoGInteraction.jl")


include("energy/energy.jl")
include("energy/energy_short_naive.jl")
include("energy/energy_long_naive.jl")
include("energy/energy_long_mid.jl")
include("energy/energy_long_wide.jl")

end

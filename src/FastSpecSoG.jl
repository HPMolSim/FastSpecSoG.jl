module FastSpecSoG

using ExTinyMD, LinearAlgebra, SpecialFunctions, ChebParticleMesh

export USeriesPara, U_series, BSA
export FSSoG_naive
export short_energy_naive, long_energy_naive, energy_naive
export mid_paras_gen, energy_mid
export FFCT_precompute, inverse_mat
export interpolate_nu_single!, interpolate_nu!

include("types.jl")

include("U_series.jl")
include("FSSoGInteraction.jl")


include("energy/energy.jl")
include("energy/energy_short_naive.jl")
include("energy/energy_long_naive.jl")
include("energy/energy_mid.jl")
include("energy/energy_long.jl")

include("FFCT/precompute.jl")
include("FFCT/nugird_ops.jl")
include("FFCT/linear_eqs.jl")

end

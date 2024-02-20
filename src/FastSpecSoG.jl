module FastSpecSoG

using ExTinyMD, LinearAlgebra, SpecialFunctions, ChebParticleMesh, SumOfExpVPMR, LoopVectorization, OMEinsum

export USeriesPara, U_series, BSA
export FSSoG_naive
export short_energy_naive, long_energy_naive, energy_naive
export mid_paras_gen, energy_mid
export FFCT_precompute, boundaries!, inverse_mat, interpolate_nu_single!, interpolate_nu!, real2Cheb!, solve_eqs!, gather_nu, gather_nu_single, zeroth_order
export energy_long, long_energy_us, long_energy_sw_0, long_energy_sw_k
export SoePara4, SoePara8, SoePara16

include("types.jl")

include("U_series.jl")
include("FSSoGInteraction.jl")


include("energy/energy.jl")
include("energy/energy_short_naive.jl")
include("energy/energy_long_naive.jl")
include("energy/energy_mid.jl")
include("energy/energy_long.jl")

include("FFCT/precompute.jl")
include("FFCT/nugrid_interpolate.jl")
include("FFCT/linear_eqs.jl")
include("FFCT/nugrid_gather.jl")
include("FFCT/zeroth_order.jl")

end

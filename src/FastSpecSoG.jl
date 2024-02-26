module FastSpecSoG

using ExTinyMD, LinearAlgebra, SpecialFunctions, ChebParticleMesh, SumOfExpVPMR, LoopVectorization, OMEinsum

export USeriesPara, U_series, BSA
export FSSoG_naive
export short_energy_naive, long_energy_naive, energy_naive
export mid_paras_gen, energy_mid
export FFCT_precompute, boundaries!, inverse_mat
export interpolate_nu_loop_single!, interpolate_nu_loop!, interpolate_nu_einsum!, interpolate_nu_einsum_non_inplace!
export real2Cheb!, solve_eqs!
export ather_nu, gather_nu_single
export zeroth_order
export energy_long_loop_k, energy_long_einsum_k
export long_energy_us_k, long_energy_us_0, long_energy_sw_0, long_energy_sw_k
export energy_long_0
export SoePara, SoePara4, SoePara8, SoePara16

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

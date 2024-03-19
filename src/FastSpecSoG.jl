module FastSpecSoG

using ExTinyMD, LinearAlgebra, SpecialFunctions, ChebParticleMesh, SumOfExpVPMR, LoopVectorization, OMEinsum, FastChebInterp

import FastChebInterp: ChebPoly

export USeriesPara, U_series, BSA, ChebUSeries
export FSSoG_naive, FSSoGInteraction
export short_energy_naive, long_energy_naive, energy_naive

export Es_Cheb_precompute, short_energy_Cheb

export mid_paras_gen, energy_mid, energy_long, energy_short
export FFCT_precompute, boundaries!, inverse_mat
export interpolate_nu_einsum!, interpolate_nu_cheb!, interpolate_nu_loop!
export real2Cheb!, solve_eqs!, gather_nu

export energy_long_loop_k, energy_long_einsum_k, energy_long_cheb_k, energy_long_0

export long_energy_us_k, long_energy_us_0, long_energy_sw_0, long_energy_sw_k
export SoePara, SoePara4, SoePara8, SoePara16

include("types.jl")

include("U_series.jl")
include("FSSoGInteraction.jl")

include("FFCT/precompute.jl")
include("FFCT/nugrid_interpolate.jl")
include("FFCT/linear_eqs.jl")
include("FFCT/nugrid_gather.jl")
include("FFCT/zeroth_order.jl")
include("FFCT/cheb_interpolate.jl")

include("energy/energy_short.jl")
include("energy/energy_mid.jl")
include("energy/energy_long.jl")

include("energy/energy.jl")
include("energy/energy_short_naive.jl")
include("energy/energy_long_naive.jl")

end

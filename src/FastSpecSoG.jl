module FastSpecSoG

using ExTinyMD, LinearAlgebra, SpecialFunctions, ChebParticleMesh, SumOfExpVPMR, LoopVectorization, OMEinsum, FastChebInterp, Polynomials, FFTW

import FastChebInterp: ChebPoly

export USeriesPara, U_series, BSA, ChebUSeries, proper_M
export FSSoG_naive, FSSoGInteraction, FSSoGThinInteraction
export short_energy_naive, long_energy_naive, energy_naive

export Es_Cheb_precompute, short_energy_Cheb

export mid_paras_gen, energy_mid, energy_long, energy_short
export long_paras_gen, interpolate_nu_cheb!, interpolate_nu_loop!
export real2Cheb!, Cheb2real!, gather_nu
export energy_long_loop_k, energy_long_cheb_k, energy_long_0

export long_energy_us_k, long_energy_us_0, long_energy_sw_0, long_energy_sw_k
export SoePara, SoePara4, SoePara8, SoePara16

export TaylorUSeries_k, TaylorUSeries
export thin_paras_gen, interpolate_thin!, gather_thin, scale_thin!
export piecewise_fft!, piecewise_ifft!, piecewise_mul!
export energy_long_thin_k

include("types.jl")

include("U_series.jl")
include("FSSoGInteraction.jl")

include("FFCT/precompute.jl")
include("FFCT/interpolate.jl")
include("FFCT/gather.jl")
include("FFCT/zeroth_order.jl")
include("FFCT/Chebyshev.jl")
include("FFCT/Taylor.jl")
include("FFCT/scaling.jl")

include("energy/energy_short.jl")
include("energy/energy_mid.jl")
include("energy/energy_long.jl")

include("energy/energy.jl")
include("energy/energy_short_naive.jl")
include("energy/energy_long_naive.jl")

end

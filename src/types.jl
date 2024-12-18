struct USeriesPara{T}
    sw::Vector{NTuple{2, T}}
end

Base.show(io::IO, uspara::USeriesPara{T}) where{T} = print(io, "USeriesPara($T), M=$(length(uspara.sw))")

struct FSSoG_naive{T} <: ExTinyMD.AbstractInteraction
    b::T
    σ::T
    ω::T
    M::Int
    r_c::T
    k_c::T
    ϵ::T
    L::NTuple{3, T}
    k_set::Vector{NTuple{3, T}}
    uspara::USeriesPara{T}
    n_atoms::Int
end

Base.show(io::IO, fssog_naive::FSSoG_naive{T}) where{T} = print(io, "FSSoG_naive($T), b=$(fssog_naive.b), σ=$(fssog_naive.σ), ω=$(fssog_naive.ω), M=$(fssog_naive.M), r_c=$(fssog_naive.r_c), k_c=$(fssog_naive.k_c), ϵ=$(fssog_naive.ϵ), L=$(fssog_naive.L), n_atoms=$(fssog_naive.n_atoms)")

mutable struct FSSoGInteraction{T} <: ExTinyMD.AbstractInteraction
    b::T
    σ::T
    ω::T
    M::Int
    ϵ::T
    L::NTuple{3, T}
    boundary::Boundary{T}
    n_atoms::Int

    uspara::USeriesPara{T}
    position::Vector{NTuple{3, T}}
    charge::Vector{T}

    # for short range
    uspara_cheb::ChebPoly{1, T, T}
    r_c::T
    F0::T

    # for 3D FFT
    gridinfo::GridInfo{3, T}
    gridbox::GridBox{3, T}
    cheb_coefs::NTuple{3, ChebCoef{T}}
    scalefactor::ScalingFactor{3, T}

    # for FFCT (using the energy_long_cheb as defalut)
    M_mid::Int
    k_x::Vector{T}
    k_y::Vector{T}
    r_z::Vector{T}
    phase_x::Vector{Complex{T}}
    phase_y::Vector{Complex{T}}
    H_r::Array{Complex{T}, 3}
    H_c::Array{Complex{T}, 3}
    cheb_mat::Array{ChebPoly{1, T, T}, 2}

    # for FFCT zero order term
    Q_0::Int
    r_z0::Vector{T}
    chebcoefs0::Vector{T}
    grids0::Vector{T}
    chebuseries::ChebPoly{1, T, T}
end

Base.show(io::IO, fssog::FSSoGInteraction{T}) where{T} = print(io, "FSSoGInteraction($T), b=$(fssog.b), σ=$(fssog.σ), ω=$(fssog.ω), M=$(fssog.M), M_mid=$(fssog.M_mid), ϵ=$(fssog.ϵ), L=$(fssog.L), r_c=$(fssog.r_c), Q_0=$(fssog.Q_0), N_grid_mid=$(fssog.gridinfo.N_pad), N_grid_long=$(size(fssog.H_r))")

mutable struct FSSoGThinInteraction{T} <: ExTinyMD.AbstractInteraction
    b::T
    σ::T
    ω::T
    M::Int
    ϵ::T
    L::NTuple{3, T}
    boundary::Boundary{T}
    n_atoms::Int

    uspara::USeriesPara{T}
    position::Vector{NTuple{3, T}}
    charge::Vector{T}

    # for short range
    uspara_cheb::ChebPoly{1, T, T}
    r_c::T
    F0::T

    # for long range by 2dfft
    r_z::Vector{T}
    H_r::Array{Complex{T}, 3}
    H_c::Array{Complex{T}, 3}
    gridinfo::GridInfo{2, T}
    pad_grids::Vector{Array{Complex{T}, 3}}
    scalefactors::Vector{ScalingFactor{2, T}}
    cheb_coefs::NTuple{2, ChebCoef{T}}
    cheb_value::Vector{Array{T, 1}}

    # for zero order term
    Q_0::Int
    r_z0::Vector{T}
    chebcoefs0::Vector{T}
    grids0::Vector{T}
    chebuseries::ChebPoly{1, T, T}
end

Base.show(io::IO, fssog_thin::FSSoGThinInteraction{T}) where{T} = print(io, "FSSoGThinInteraction($T), b=$(fssog_thin.b), σ=$(fssog_thin.σ), ω=$(fssog_thin.ω), M=$(fssog_thin.M), ϵ=$(fssog_thin.ϵ), L=$(fssog_thin.L), r_c=$(fssog_thin.r_c), Q_0=$(fssog_thin.Q_0), N_grid=$(size(fssog_thin.H_r))")

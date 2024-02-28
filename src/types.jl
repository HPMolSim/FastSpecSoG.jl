struct USeriesPara{T}
    sw::Vector{NTuple{2, T}}
end

struct TaylorSeriesPara{N, T}
    coefs::NTuple{N, T}
end

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

mutable struct FSSoGInteraction{T} <: ExTinyMD.AbstractInteraction
    b::T
    σ::T
    ω::T
    M::Int
    ϵ::T
    L::NTuple{3, T}
    uspara::USeriesPara{T}
    n_atoms::Int

    # para for 3D FFT
    n_k::NTuple{3, Int}
    h::NTuple{3, T}
    win_cut::Int
    win_width::NTuple{3, T} # (win_cut + 0.5) * h
    β::T # 5 * win_cut

    pad_ratio::Int
    pad_cut::Int # pad_ratio * win_cut
    pad_width::T # pad_cut * h

    kx::Vector{T}
    ky::Vector{T}
    kz::Vector{T}
    TdK::Array{T, 3}
    grid_extended::Array{T, 3}
    grid_real::Array{T, 3}
    Pkbx::Array{T, 2}
    Pkby::Array{T, 2}
    Pkbz::Array{T, 2}


end
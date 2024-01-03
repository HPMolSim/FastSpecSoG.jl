struct USeriesPara{T}
    sw::Vector{NTuple{2, T}}
end

const preset_parameters = [
    [2.0, 5.027010924194599, 0.994446492762232252, 6], 
    [1.62976708826776469, 3.633717409009413, 1.00780697934380681, 16], 
    [1.48783512395703226, 2.662784519725113, 0.991911705759818, 30], 
    [1.32070036405934420, 2.277149356440992, 1.00188914114811980, 64],
    [1.21812525709410644, 1.774456369233284, 1.00090146156033341, 102]]

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
struct FastSpecSOG{T} <: ExTinyMD.AbstractInteraction
    b::T
    σ::T
    M::Int
    r_c::T
end

struct USeriesPara{T}
    sw::Vector{NTuple{2, T}}
end
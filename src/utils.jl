function rmsd(y1::Vector{T}, y2::Vector{T}) where T
    @assert length(y1) == length(y2)
    return sqrt(sum((y1 - y2).^2) / length(y1))
end

function rrmsd(y1::Vector{T}, y2::Vector{T}) where T
    @assert length(y1) == length(y2)
    return rmsd(y1, y2) / abs(sum(y2) / length(y2))
end
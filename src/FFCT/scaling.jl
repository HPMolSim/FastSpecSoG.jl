function scale_thin!(pad_grids::Vector{Array{Complex{T}, 3}}, scalefactors::Vector{ScalingFactor{2, T}}, H_r::Array{Complex{T}, 3}) where{T}
    for i in 1:length(pad_grids)
        piecewise_mul!(pad_grids[i], scalefactors[i].factors)
    end
    return sum!(H_r, pad_grids)
end

function piecewise_mul!(grid::Array{Complex{T}, 3}, factor::Array{Complex{T}, 2}) where{T}
    for k in 1:size(grid, 3)
        (@view grid[:, :, k]) .*= factor
    end
    return grid
end

function piecewise_fft!(grid::Array{Complex{T}, 3}) where{T}
    for k in 1:size(grid, 3)
        fft!(@view grid[:, :, k])
    end
    return grid
end

function piecewise_ifft!(grid::Array{Complex{T}, 3}) where{T}
    for k in 1:size(grid, 3)
        ifft!(@view grid[:, :, k])
    end
    return grid
end

function sum!(pad_gird0::Array{Complex{T}, 3}, pad_grids::Vector{Array{Complex{T}, 3}}) where{T}
    set_zeros!(pad_gird0)
    for i in 1:length(pad_grids)
        pad_gird0 .+= pad_grids[i]
    end
    return pad_gird0
end
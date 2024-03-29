function scale_thin!(gridboxs::Array{GridBox{2, T}, 2}, scalefactors::Vector{ScalingFactor{2, T}}, H_r::Array{Complex{T}, 3}) where{T}
    set_zeros!(H_r)
    for i in 1:size(gridboxs, 1)
        for j in 1:size(gridboxs, 2)
            t = fft(gridboxs[i, j].pad_grid)
            (@view H_r[:, :, i]) .+= t .* scalefactors[j].factors
        end
    end
    return H_r
end

function scaling_Q2D!(grid::Array{Complex{T}, 3}, factor::Array{Complex{T}, 2}) where{T}
    for k in 1:size(grid, 3)
        (@view grid[:, :, k]) .*= factor
    end
    return grid
end

function fft_Q2d!(grid::Array{Complex{T}, 3}) where{T}
    for k in 1:size(grid, 3)
        fft!(@view grid[:, :, k])
    end
    return grid
end

function ifft_Q2d!(grid::Array{Complex{T}, 3}) where{T}
    for k in 1:size(grid, 3)
        ifft!(@view grid[:, :, k])
    end
    return grid
end
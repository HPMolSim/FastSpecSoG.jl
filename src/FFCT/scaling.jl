function thin_scale!(gridboxs::Array{GridBox{2, T}, 2}, scalefactors::Vector{ScalingFactor{T}}) where{T}
    for i in 1:size(gridboxs, 1)
        for j in 1:size(gridboxs, 2)
            scale!(fft!(gridboxs[i, j]), scalefactors[j])
        end
    end
    return gridboxs
end

function sum_gridboxs!(gridboxs::Array{GridBox{2, T}, 2}, H_r::Array{Complex{T}, 3}) where{T}
    set_zeros!(H_r)
    for i in 1:size(gridboxs, 1)
        for j in 1:size(gridboxs, 2)
            (@view H_r[:, :, i]) .+= gridboxs[i, j].pad_grid
        end
    end
    return H_r
end

function ifft2d!(H_c::Array{Complex{T}, 3}) where{T}
    for i in 1:size(H_c, 3)
        ifft!(@view(H_c[:, :, i]))
    end
    return H_c
end
using FINUFFT

function chebyshev_2d_grid(Lx::T, Ly::T, Nx::Int, Ny::Int) where{T}
    rx = Vector{T}([Lx / 2 * ( 1.0 - cos((2i - 1)*π / (2 * Nx)) ) for i in 1:Nx])
    ry = Vector{T}([Ly / 2 * ( 1.0 - cos((2i - 1)*π / (2 * Ny)) ) for i in 1:Ny])
    return rx, ry
end

function eval_2d_chebgrid!(chebgrid::AbstractArray{T, 2}, f::Function, rx::Vector{T}, ry::Vector{T}) where{T}
    @inbounds for i in 1:length(rx)
        for j in 1:length(ry)
            chebgrid[i, j] = f(rx[i], ry[j])
        end
    end
    return chebgrid
end

eval_2d_chebgrid(f::Function, rx::Vector{T}, ry::Vector{T}) where{T} = eval_2d_chebgrid!(zeros(T, length(rx), length(ry)), f, rx, ry)
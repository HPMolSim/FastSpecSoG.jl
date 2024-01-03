# in this script we will define a few functions to be used

@inline function horner(x::T, p::AbstractArray{T, 1}) where{T}
    ex = p[end]
    for i = length(p)-1:-1:1
        ex = p[i] + x * ex
    end
    return ex
end

function ChebPoly(x::T, scale::T, n::TI) where{T<:Real, TI<:Integer}
    return T(cos(T(n) * acos(x / scale)))
end

"""
function Wkb(x::T, width::T, β::T) where{T<:Real}

    Wkb kernel function
"""
function Wkb(x::T, width::T, β::T) where{T<:Real}
    return T(besseli(0, β * sqrt(one(T) - (x / width)^2)) / besseli(0, β) * (abs(x) <= width))
end

function Wkb(x::Vector{T}, width::T, β::T) where{T<:Real}
    return T.(besseli.(0, β .* sqrt.(one(T) .- (x ./ width).^2)) ./ besseli(0, β) .* (abs.(x) .<= width))
end

"""
function Pkb_coef(h::T, wsearch::TI; P::TI = TI(2 * wsearch), width::T = T(wsearch * h), β::T = T(wsearch * 5), ν::TI = 10) where{T <: Real, TI <: Integer}

    h: grid size
    wsearch: the number of grid points in the search window
    P: window number, default: 2 * wsearch
    width: the width of the window, default: wsearch * h
    β: the parameter in the WKB function, default: 5 * wsearch
    ν: the number of Chebyshev polynomials, default: 10

this function will return the coefficients for the WKB function in form of parameters of the Chebyshev polynomials in each interval
"""
function Pkb_coef(h::T, wsearch::TI; P::TI = TI(2 * wsearch), width::T = T(wsearch * h), β::T = T(wsearch * 5), ν::TI = 10) where{T <: Real, TI <: Integer}
    width_window = T(width + h / 2)
    C = zeros(T, P + 2, ν)

    x = chebyshevpoints(T, ν, Val(1))

    for i in 1:P + 2
        if i == 1
            center = - width - h / 4
            scale = h / 4
        elseif i == P + 2
            center = width + h / 4
            scale = h / 4
        else
            center = - width + (i - 1.5) * h
            scale = h / 2
        end

        b = Wkb(center .+ scale .* x, width_window, β)

        C[i, :] = chebyshevtransform(b, Val(1))
    end

    return C
end

function Pkb_cal(x::T, C::AbstractArray{T, 1}, i::TI, h::T, wsearch::TI; width::T = T(wsearch * h)) where{T <: Real, TI <: Integer}
    if i == 1
        center = - width - h / 4
        scale = h / 4
    elseif i == length(C)
        center = width + h / 4
        scale = h / 4
    else
        center = - width + (i - 1.5) * h
        scale = h / 2
    end

    x = (x - center) / scale

    return horner(x, C)
end

function FWkb(k::T, width::T, β::T) where{T}
    return T(2 * width * sinh(sqrt(β^2 - (k * width)^2)) / (besseli(0, β) * sqrt(β^2 - (k * width)^2)))
end
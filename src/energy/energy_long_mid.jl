function Tdk(kx::T, ky::T, kz::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    k2 = kx^2 + ky^2 + kz^2

    Tdk_value = zero(T)
    for i in 1:M_mid
        sl, wl = uspara[i]
        Tdk_value += wl * sl^3 * exp(- sl^2 * k2 / 4)
    end

    return sqrt(π)/4 * Tdk_value
end

function mid_paras_gen(N_real::NTuple{3, Int}, w::NTuple{3, Int}, L::NTuple{3, T}, extra_pad_ratio::Int, cheb_order::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    periodicity = (true, true, false)
    extra_pad = extra_pad_ratio .* w
    @assert M_mid <= length(uspara)

    gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
    gridbox = GridBox(gridinfo)

    f_window = [x -> Wkb(x, (w[i] + 0.5) * gridinfo.h[i], 5.0 * w[i]) for i in 1:3]
    cheb_coefs = tuple([ChebCoef(f_window[i], gridinfo.h[i], w[i], cheb_order) for i in 1:3]...)

    F_f_window = [x -> FWkb(x, (w[i] + 0.5) * gridinfo.h[i], 5.0 * w[i]) for i in 1:3]

    func_scale = (kx, ky, kz) -> (F_f_window[1](kx) * F_f_window[2](ky) * F_f_window[3](kz))^(-2) * Tdk(kx, ky, kz, uspara, M_mid)
    scalefactor = ScalingFactor(func_scale, gridinfo)

    return (gridinfo, gridbox, cheb_coefs, scalefactor)
end

function energy_long_mid(qs::Vector{T}, poses::Vector{NTuple{3, T}}, gridinfo::GridInfo{3, T}, gridbox::GridBox{3, T}, cheb_coefs::NTuple{N, ChebCoef{T}}, scalefactor::ScalingFactor{3, T}) where{T}

    interpolate!(qs, poses, gridinfo, gridbox, cheb_coefs)
    fft!(gridbox)
    scale!(gridbox, scalefactor)
    ifft!(gridbox)
    E = gather(qs, poses, gridinfo, gridbox, cheb_coefs)

    return E
end
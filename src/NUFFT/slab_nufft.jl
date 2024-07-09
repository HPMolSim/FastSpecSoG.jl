using FINUFFT
export SlabNUFFT, nufft_energy_long_k

function Tdk_M(kx::T, ky::T, uspara::USeriesPara{T}, M::Int) where{T}
    k2 = kx^2 + ky^2
    if iszero(k2)
        return zero(T)
    else
        sl, wl = uspara.sw[M]
        Tdk_value = wl * sl^2 * exp(- sl^2 * k2 / 4)
        return Tdk_value
    end
end

@inbounds function real2Cheb_nufft!(f_real::Array{Complex{T}, 3}, f_cheb::Array{Complex{T}, 3}, r_z::Vector{T}, L_z::T) where{T}

    set_zeros!(f_cheb)
    N_cheb = size(f_real, 1)

    for k in 1:N_cheb, l in 1:N_cheb
        cheb_temp = k == 1 ? 0.5 : chebpoly(k - 1, r_z[l] - L_z / T(2), L_z / T(2))
        for j in 1:size(f_real, 3), i in 1:size(f_real, 2)
            f_cheb[k, i, j] += 2 / N_cheb * f_real[l, i, j] * cheb_temp
        end
    end

    return f_cheb
end

mutable struct SlabNUFFT{T}
    N_k::Tuple{Int, Int}
    N_cheb::Int
    L::Tuple{T, T, T}
    eps::T
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    q::Vector{Complex{T}}
    rz::Vector{T}

    uspara::USeriesPara{T}
    M_mid::Int

    scaling_factor::Vector{Array{T, 2}} # for the M-th order of sog
    f_real::Array{Complex{T}, 3}
    f_cheb::Array{Complex{T}, 3}
    f_temp::Array{Complex{T}, 2}
    ϕ::Array{Complex{T}, 2}
    ϕ_temp::Array{Complex{T}, 1}
end

function SlabNUFFT(n_atoms::Int, N_k::Tuple{Int, Int}, N_cheb::Int, L::Tuple{T, T, T}, eps::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    x = zeros(T, n_atoms)
    y = zeros(T, n_atoms)
    z = zeros(T, n_atoms)
    q = zeros(Complex{T}, n_atoms)
    rz = Vector{T}([L[3] / 2 * ( 1.0 - cos((2i - 1)*π / (2 * N_cheb)) ) for i in 1:N_cheb])

    k_x = [2π * i / L[1] for i in Int(ceil(-N_k[1]/2)):Int(floor((N_k[1] - 1)/2))]
    k_y = [2π * i / L[2] for i in Int(ceil(-N_k[2]/2)):Int(floor((N_k[2] - 1)/2))]
    scaling_factor = [zeros(T, N_k...) for _ in M_mid + 1:length(uspara.sw)]
    for m in M_mid + 1:length(uspara.sw)
        for i in 1:N_k[1]
            for j in 1:N_k[2]
                scaling_factor[m - M_mid][i, j] = Tdk_M(k_x[i], k_y[j], uspara, m)
            end
        end
    end

    # f = [zeros(Complex{T}, N_k...) for i in M_mid + 1:length(uspara.sw), j in 1:N_cheb]
    f_real = zeros(Complex{T}, N_cheb, N_k[1], N_k[2])
    f_cheb = zeros(Complex{T}, N_cheb, N_k[1], N_k[2])
    f_temp = zeros(Complex{T}, N_k...)
    ϕ = zeros(Complex{T}, N_cheb, n_atoms)
    ϕ_temp = zeros(Complex{T}, n_atoms)

    return SlabNUFFT(N_k, N_cheb, L, eps, x, y, z, q, rz, uspara, M_mid, scaling_factor, f_real, f_cheb, f_temp, ϕ, ϕ_temp)
end

function update_xyz_ffct!(qs::Vector{T}, poses::Vector{NTuple{3, T}}, slab_nufft::SlabNUFFT{T}) where{T}
    for i in 1:length(qs)
        slab_nufft.x[i] = 2π * poses[i][1] / slab_nufft.L[1]
        slab_nufft.y[i] = 2π * poses[i][2] / slab_nufft.L[2]
        slab_nufft.z[i] = poses[i][3]
    end
    return slab_nufft
end

function update_q_ffct!(qs::Vector{T}, M::Int, rz::T, slab_nufft::SlabNUFFT{T}) where{T}
    n_atoms = length(qs)
    for i in 1:n_atoms
        z = slab_nufft.z[i]
        sl, wl = slab_nufft.uspara.sw[M]
        qrzM = qs[i] * exp(-(z - rz)^2 / sl^2)
        slab_nufft.q[i] = Complex{T}(qrzM)
    end
    return slab_nufft
end

function nufft_energy_long_k(qs::Vector{T}, poses::Vector{NTuple{3, T}}, slab_nufft::SlabNUFFT{T}) where{T}

    n_atoms = length(qs)
    update_xyz_ffct!(qs, poses, slab_nufft)
    set_zeros!(slab_nufft.f_real)
    for n_cheb in 1:slab_nufft.N_cheb
        for l in slab_nufft.M_mid + 1:length(slab_nufft.uspara.sw)
            update_q_ffct!(qs, l, slab_nufft.rz[n_cheb], slab_nufft)
            nufft2d1!(slab_nufft.x, slab_nufft.y, slab_nufft.q, -1, slab_nufft.eps, slab_nufft.f_temp)
            (@view slab_nufft.f_real[n_cheb, :, :]) .+= slab_nufft.f_temp .* slab_nufft.scaling_factor[l - slab_nufft.M_mid]
        end
    end

    real2Cheb_nufft!(slab_nufft.f_real, slab_nufft.f_cheb, slab_nufft.rz, slab_nufft.L[3])

    ϕ_temp = slab_nufft.ϕ_temp
    for n_cheb in 1:slab_nufft.N_cheb
        slab_nufft.f_temp .= @view slab_nufft.f_cheb[n_cheb, :, :]
        nufft2d2!(slab_nufft.x, slab_nufft.y, ϕ_temp, 1, slab_nufft.eps, slab_nufft.f_temp)
        (@view slab_nufft.ϕ[n_cheb, :]) .= ϕ_temp
    end

    L_z = slab_nufft.L[3]
    E = zero(Complex{T})
    for i in 1:n_atoms
        cheb_coefs = @view slab_nufft.ϕ[:, i]
        z = slab_nufft.z[i]
        for k in 1:slab_nufft.N_cheb
            cheb_val = chebpoly(k - 1, z - L_z / T(2), L_z / T(2))
            E += qs[i] * cheb_coefs[k] * cheb_val
        end
    end

    return real(E) / (slab_nufft.L[1] * slab_nufft.L[2]) * π / 2
end
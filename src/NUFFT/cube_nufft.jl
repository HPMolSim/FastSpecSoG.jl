using FINUFFT
export CubeNUFFT, nufft_energy_mid

mutable struct CubeNUFFT{T}
    N_k::Tuple{Int, Int, Int}
    L::Tuple{T, T, T}
    Lpz::T
    eps::T
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    q::Vector{Complex{T}}
    scaling_factor::Array{T, 3}
    uspara::USeriesPara{T}
    M_mid::Int

    f::Array{Complex{T}, 3}
    ϕ::Vector{Complex{T}}
end

function CubeNUFFT(n_atoms::Int, N_k::Tuple{Int, Int, Int}, L::Tuple{T, T, T}, λ::T, eps::T, uspara::USeriesPara{T}, M_mid::Int) where{T}
    x = zeros(T, n_atoms)
    y = zeros(T, n_atoms)
    z = zeros(T, n_atoms)
    q = zeros(Complex{T}, n_atoms)

    Lpz = L[3] * λ

    k_x = [2π * i / L[1] for i in Int(ceil(-N_k[1]/2)):Int(floor((N_k[1] - 1)/2))]
    k_y = [2π * i / L[2] for i in Int(ceil(-N_k[2]/2)):Int(floor((N_k[2] - 1)/2))]
    k_z = [2π * i / Lpz for i in Int(ceil(-N_k[3]/2)):Int(floor((N_k[3] - 1)/2))]
    scaling_factor = zeros(T, N_k...)
    for i in 1:N_k[1]
        for j in 1:N_k[2]
            for k in 1:N_k[3]
                scaling_factor[i, j, k] = Tdk(k_x[i], k_y[j], k_z[k], uspara, M_mid)
            end
        end
    end

    f = zeros(Complex{T}, N_k...)
    ϕ = zeros(Complex{T}, n_atoms)

    return CubeNUFFT(N_k, L, Lpz, eps, x, y, z, q, scaling_factor, uspara, M_mid, f, ϕ)
end

function update_q!(qs::Vector{T}, cube_nufft::CubeNUFFT{T}) where{T}
    cube_nufft.q = Complex{T}.(qs)
    return cube_nufft
end

function update_xyz!(poses::Vector{NTuple{3, T}}, cube_nufft::CubeNUFFT{T}) where{T}
    for i in 1:length(poses)
        cube_nufft.x[i] = 2π * poses[i][1] / cube_nufft.L[1]
        cube_nufft.y[i] = 2π * poses[i][2] / cube_nufft.L[2]
        cube_nufft.z[i] = 2π * poses[i][3] / cube_nufft.Lpz
    end
    return cube_nufft
end

function nufft_energy_mid(qs::Vector{T}, poses::Vector{NTuple{3, T}}, cube_nufft::CubeNUFFT{T}) where{T}
    update_q!(qs, cube_nufft)
    update_xyz!(poses, cube_nufft)
    nufft3d1!(cube_nufft.x, cube_nufft.y, cube_nufft.z, cube_nufft.q, -1, cube_nufft.eps, cube_nufft.f)
    cube_nufft.f .*= cube_nufft.scaling_factor
    nufft3d2!(cube_nufft.x, cube_nufft.y, cube_nufft.z,  cube_nufft.ϕ, 1, cube_nufft.eps, cube_nufft.f)
    E = dot(qs, cube_nufft.ϕ) / 2
    return real(E) / (cube_nufft.L[1] * cube_nufft.L[2] * cube_nufft.Lpz) * 4π
end
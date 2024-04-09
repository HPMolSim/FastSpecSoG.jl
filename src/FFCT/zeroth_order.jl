function proper_N(accuracy::T, η::T) where{T}
    N = 0
    t = one(T)
    if η ≥ one(T)
        return N
    end

    while t > η * accuracy
        N += 1
        t *= η / N
    end
    return N
end

function FGT1d(qs::Vector{T}, poses::Vector{NTuple{3, T}}, s::T, L_z::T, r_z::Vector{T}, chebcoefs::Vector{T}, grids::Vector{T}) where{T}
    n_atoms = length(qs)
    Q_0 = length(r_z)
    # grids = zeros(T, Q_0)
    # r_z = [L_z / 2 * ( 1.0 - cos((2i - 1)*π / (2 * Q_0)) ) for i in 1:Q_0]

    set_zeros!(grids)
    set_zeros!(chebcoefs)

    # (L_z/s)^2N / factorials(N) < 1e-16
    N = proper_N(1e-16, (L_z / s)^2)

    for i in 1:n_atoms
        q = qs[i]
        for j in 1:Q_0
            z = poses[i][3]
            t = - (z - r_z[j])^2 / s^2
            if N == 0
                grids[j] += q * exp(t)
            else
                for k in 1:N
                    grids[j] += q * special_powern(t, k)
                end
            end
        end
    end

    # chebcoefs = zeros(T, Q_0)

    for k in 1:Q_0, l in 1:Q_0
        cheb_temp = k == 1 ? 0.5 : chebpoly(k - 1, r_z[l] - L_z / T(2), L_z / T(2))
        chebcoefs[k] += 2 / Q_0 * grids[l] * cheb_temp
    end

    E_k0 = zero(T)
    for i in 1:n_atoms
        q = qs[i]
        z = poses[i][3]
        for k in 1:Q_0
            E_k0 += q * chebcoefs[k] * chebpoly(k - 1, z - L_z / T(2), L_z / T(2))
        end
    end

    return E_k0
end

function FGT1d_naive(qs::Vector{T}, poses::Vector{NTuple{3, T}}, s::T) where{T}
    n_atoms = length(qs)
    E_k0 = big(zero(T))

    for i in 1:n_atoms, j in 1:n_atoms
        qi = qs[i]
        xi, yi, zi = poses[i]
        qj = qs[j]
        xj, yj, zj = poses[j]
        E_k0 += qi * qj * exp( - big(zi - zj)^2 / s^2)
    end

    return Float64(E_k0)
end


function zeroth_order(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, uspara::USeriesPara{T}, M_mid::Int, r_z::Vector{T}, chebcoefs::Vector{T}, grids::Vector{T}) where{T}

    L_x, L_y, L_z = L

    E_k0 = zero(T)

    for l in M_mid + 1:length(uspara.sw)
        sl, wl = uspara.sw[l]
        E_k0 += wl * sl^2 * FGT1d(qs, poses, sl, L_z, r_z, chebcoefs, grids)
    end

    return π * E_k0 / (2 * L_x * L_y)
end
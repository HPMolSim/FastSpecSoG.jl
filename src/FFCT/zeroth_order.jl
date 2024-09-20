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

function FGT1d(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L_z::T, chebuseries::ChebPoly{1, T, T}, r_z::Vector{T}, chebcoefs::Vector{T}, grids::Vector{T}) where{T}
    n_atoms = length(qs)
    Q_0 = length(r_z)

    set_zeros!(grids)
    set_zeros!(chebcoefs)


    for i in 1:n_atoms
        q = qs[i]
        for j in 1:Q_0
            z = poses[i][3]
            grids[j] += q * chebuseries(abs(z - r_z[j]))
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

function FGT1d_per_atom(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L_z::T, chebuseries::ChebPoly{1, T, T}, r_z::Vector{T}, chebcoefs::Vector{T}, grids::Vector{T}) where{T}
    n_atoms = length(qs)
    Q_0 = length(r_z)

    set_zeros!(grids)
    set_zeros!(chebcoefs)


    for i in 1:n_atoms
        q = qs[i]
        for j in 1:Q_0
            z = poses[i][3]
            grids[j] += q * chebuseries(abs(z - r_z[j]))
        end
    end

    # chebcoefs = zeros(T, Q_0)

    for k in 1:Q_0, l in 1:Q_0
        cheb_temp = k == 1 ? 0.5 : chebpoly(k - 1, r_z[l] - L_z / T(2), L_z / T(2))
        chebcoefs[k] += 2 / Q_0 * grids[l] * cheb_temp
    end

    E_k0_N = zeros(T, n_atoms)
    for i in 1:n_atoms
        q = qs[i]
        z = poses[i][3]
        for k in 1:Q_0
            E_k0_N[i] += q * chebcoefs[k] * chebpoly(k - 1, z - L_z / T(2), L_z / T(2))
        end
    end

    return E_k0_N
end


function zeroth_order(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, chebuseries::ChebPoly{1, T, T}, r_z::Vector{T}, chebcoefs::Vector{T}, grids::Vector{T}) where{T}

    L_x, L_y, L_z = L

    E_k0 = FGT1d(qs, poses, L_z, chebuseries, r_z, chebcoefs, grids)

    return E_k0 / (2 * L_x * L_y)
end

function zeroth_order_per_atom(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, chebuseries::ChebPoly{1, T, T}, r_z::Vector{T}, chebcoefs::Vector{T}, grids::Vector{T}) where{T}

    L_x, L_y, L_z = L

    E_k0_per_atom = FGT1d_per_atom(qs, poses, L_z, chebuseries, r_z, chebcoefs, grids)

    return E_k0_per_atom ./ (2 * L_x * L_y)
end
# Shared helpers for the test suite. Included first by test/runtests.jl, before
# any test file that uses them.

# An array-of-structs element that is NOT an NTuple and NOT an SVector, and
# supports nothing but indexing. If any read-only position parameter has
# crept back to a concrete `Vector{NTuple{3,T}}` annotation, or if a kernel
# destructures a position instead of indexing it, every query below fails.
struct IdxOnlyPoint
    v::NTuple{3, Float64}
end
Base.getindex(p::IdxOnlyPoint, i::Int) = p.v[i]

# Squared slab minimum-image distance, written out independently of
# `FastSpecSoG._min_image_slab` so that the configurations below are not built
# with the same code they are used to test.
function _slab_r2(p, q, L)
    dx = p[1] - q[1]; dx -= L[1] * round(dx / L[1])
    dy = p[2] - q[2]; dy -= L[2] * round(dy / L[2])
    dz = p[3] - q[3]
    return dx^2 + dy^2 + dz^2
end

"""
    _random_poses(n, L; min_r, initial)

`n` further positions uniform in the box, each at least `min_r` from every
position already placed (including those in `initial`). A minimum separation is
required, not cosmetic: the short-range Chebyshev interpolant is only built on
`[r_min, r_c]`, so a pair closer than `r_min` raises an out-of-domain error from
FastChebInterp. This is the framework-free equivalent of the `min_r = 1.0`
that `SimulationInfo` applies in the older tests.
"""
function _random_poses(n, L; min_r = 1.0, initial = NTuple{3, Float64}[],
                       max_attempts = 10_000)
    poses = copy(initial)
    target = length(poses) + n
    while length(poses) < target
        placed = false
        for _ in 1:max_attempts
            p = (rand() * L[1], rand() * L[2], rand() * L[3])
            if all(q -> _slab_r2(p, q, L) >= min_r^2, poses)
                push!(poses, p)
                placed = true
                break
            end
        end
        placed || error("could not place particle $(length(poses) + 1) at min_r = $min_r")
    end
    return poses
end

# ---------------------------------------------------------------------------
# A per-atom quasi-2D Ewald reference
# ---------------------------------------------------------------------------
#
# ExTinyMD 0.3 exposes `coulomb_energy`, `short_energy` and `long_energy`,
# all of which return totals. The EwaldSummations functions this suite used to
# call, `Ewald2D_short_energy_N` and `Ewald2D_long_energy_N`, returned per-atom
# VECTORS, so the per-atom comparison in test/energy.jl has no direct
# replacement in the new API and is written out here.
#
# It is transcribed from the same formulas ExTinyMD implements
# (src/interactions/electrostatics/short.jl and long_ewald2d.jl) but is not
# vacuous: test/energy.jl requires its short part to equal `short_energy`, its
# long part to equal `long_energy`, and its total to equal `coulomb_energy`, so
# a mistake here shows up as a failure rather than as a tolerance that passes.
#
# Crucially, `α` and `s` are taken as keywords here and `r_c`/`k_c` are
# recomputed as `s/α` and `2αs` locally, NOT read off the interaction. So this
# reference does not follow a swapped α/s pair: if the `Ewald2D` call at the
# test site had them the wrong way round, the assertions against this reference
# would fail even though the total energy is α-independent.

const _EXP_OVERFLOW = log(floatmax(Float64)) - 1.0

# exp(kz) * erfc(arg), guarded: for large positive kz the exp overflows while
# the erfc underflows, and the product tends to zero, so Inf * 0 must not be
# allowed to produce NaN.
function _exp_erfc(kz::Float64, arg::Float64)
    kz > _EXP_OVERFLOW && return 0.0
    return exp(kz) * erfc(arg)
end

"""
    ewald2d_per_atom(poses, charges, L; α, s, ϵ = 1.0) -> (E_short, E_long)

Per-atom real-space and reciprocal-space quasi-2D Ewald energies. `O(N^2 K)`.
x and y are periodic under `L[1]`/`L[2]`; z is free.
"""
function ewald2d_per_atom(poses, charges, L; α, s, ϵ = 1.0)
    n = length(charges)
    r_c = s / α
    k_c = 2 * α * s
    A = L[1] * L[2]

    # real space: half of each pair term to each partner, plus the self term
    E_short = zeros(Float64, n)
    for i in 1:n - 1, j in i + 1:n
        r = sqrt(_slab_r2(poses[i], poses[j], L))
        (r < r_c && r > 0.0) || continue
        t = charges[i] * charges[j] * erfc(α * r) / r
        E_short[i] += t / 2
        E_short[j] += t / 2
    end
    for i in 1:n
        E_short[i] -= charges[i]^2 * α / sqrt(π)
    end
    E_short ./= (4π * ϵ)

    E_long = zeros(Float64, n)
    # k = 0 term
    for i in 1:n, j in 1:n
        z = poses[i][3] - poses[j][3]
        E_long[i] -= charges[i] * charges[j] *
                     (exp(-(α * z)^2) / (α * sqrt(π)) + z * erf(α * z)) / (4 * A)
    end
    # k != 0 terms, over the same reciprocal set ExTinyMD builds (0 < |k| <= k_c)
    mx_max = ceil(Int, k_c * L[1] / 2π) + 1
    my_max = ceil(Int, k_c * L[2] / 2π) + 1
    for mx in -mx_max:mx_max, my in -my_max:my_max
        kx = 2π * mx / L[1]
        ky = 2π * my / L[2]
        k = sqrt(kx^2 + ky^2)
        (0 < k <= k_c) || continue
        for i in 1:n
            acc = 0.0
            for j in 1:n
                x = poses[i][1] - poses[j][1]
                y = poses[i][2] - poses[j][2]
                z = poses[i][3] - poses[j][3]
                g = _exp_erfc(k * z, k / (2α) + α * z) +
                    _exp_erfc(-k * z, k / (2α) - α * z)
                acc += charges[i] * charges[j] * cos(kx * x + ky * y) * g
            end
            E_long[i] += acc / (8 * A * k)
        end
    end
    E_long ./= ϵ

    return E_short, E_long
end

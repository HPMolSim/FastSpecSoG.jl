function F0_cal(b::T, σ::T, ω::T, M::Int) where{T}
    return - log(b) / sqrt(2π) / σ * (ω + (one(T) - b^(-M)) / (b - one(T)))
end

function F0_cal(;preset::Int = 1, r_c::T = 10.0) where{T}
    @assert preset ≤ length(preset_parameters)

    b, σ, ω, M = T(preset_parameters[preset][1]), T(preset_parameters[preset][2]) * r_c, T(preset_parameters[preset][3]), Int(preset_parameters[preset][4])

    return F0_cal(b, σ, ω, M)
end

function Es_USeries_Cheb(uspara::USeriesPara{T}, r_min::T, r_max::T, Q::Int) where{T}
    f = r -> U_series(r, uspara)
    x = chebpoints(Q, r_min, r_max)

    return chebinterp(f.(x), r_min, r_max)
end

function Es_Cheb_precompute(preset::Int, r_min::T, r_max::T, Q::Int) where{T}
    uspara = USeriesPara(preset, r_c = r_max)
    uspara_cheb = Es_USeries_Cheb(uspara, r_min, r_max, Q)
    F0 = F0_cal(preset = preset, r_c = r_max)
    return uspara_cheb, F0
end

function Es_self(q::Vector{T}, F0::T) where{T}

    Q = sum(qi^2 for qi in q)

    return Q * F0
end

function Es_Cheb_pair(q_1::T, q_2::T, uspara_cheb::ChebPoly{1, T, T}, r::T) where{T}
    return q_1 * q_2 * (one(T) / r - uspara_cheb(r))
end

# ---------------------------------------------------------------------------
# Candidate pairs for the short-range sums
# ---------------------------------------------------------------------------
#
# A supplied `neighbor_list` is treated as CANDIDATE PAIRS ONLY. Only `pair[1]`
# and `pair[2]` are read; any distance the list carries is ignored and the true
# three-dimensional separation is recomputed from the positions by
# `_min_image_slab`. That is not defensive padding: a quasi-2D cell list builds
# its list over in-plane coordinates and so reports the in-plane distance, and
# an all-pairs finder reports zero. Trusting a supplied `r` is what made
# ExTinyMD's own `Ewald2D` short-range sum return +0.0238 against a true
# -0.1539 before Phase 1 fixed it.
#
# An in-plane candidate list is always a superset of the true pair set (in-plane
# distance <= 3-D distance), so re-filtering on the recomputed distance recovers
# the correct set rather than merely rejecting bad data.
#
# With no list at all, every i < j pair is a candidate: O(N^2), correct, and
# framework-free. Pass `neighbor_list` for anything larger than a smoke test.

# Accept either a bare pair list or an object that owns one under a
# `neighbor_list` property (which is how MD-framework cell lists expose theirs),
# so a caller need not reach inside their finder.
@inline _pair_list(::Nothing) = nothing
@inline _pair_list(neighbor_list) =
    hasproperty(neighbor_list, :neighbor_list) ? neighbor_list.neighbor_list : neighbor_list

_candidate_pairs(::Nothing, n_atoms::Int) = ((i, j) for i in 1:n_atoms - 1 for j in i + 1:n_atoms)
_candidate_pairs(neighbor_list, ::Int) = neighbor_list

# Kept as its own function so that each candidate-iterator type specializes.
function _short_Cheb_sum(pairs, uspara_cheb::ChebPoly{1, T, T}, r_c::T, L::NTuple{3, T}, position, q::Vector{T}) where{T}
    r_c_sq = r_c^2
    energy_short = zero(T)

    for pair in pairs
        i, j = pair[1], pair[2]
        _, _, r_sq = _min_image_slab(position[i], position[j], L)
        # Two separate rejections, where `position_check3D`'s all-zero sentinel
        # used to conflate them into one `iszero(r_sq)` guard:
        #   r_sq >= r_c_sq -- outside the cutoff, as before;
        #   iszero(r_sq)   -- a coincident pair. `Es_Cheb_pair` computes 1/r,
        #                     so keeping these would produce Inf (and 0/0 = NaN
        #                     once multiplied by a zero charge). They arise
        #                     from any lattice initialisation: two particles at
        #                     the same site, or separated by an exact multiple
        #                     of Lx or Ly with the same y and z. The old
        #                     sentinel skipped them incidentally; this skips
        #                     them deliberately.
        (r_sq < r_c_sq && !iszero(r_sq)) || continue
        energy_short += Es_Cheb_pair(q[i], q[j], uspara_cheb, sqrt(r_sq))
    end

    return energy_short
end

function _short_Cheb_sum_per_atom!(energy_short_per_atoms::Vector{T}, pairs, uspara_cheb::ChebPoly{1, T, T}, r_c::T, L::NTuple{3, T}, position, q::Vector{T}) where{T}
    r_c_sq = r_c^2

    for pair in pairs
        i, j = pair[1], pair[2]
        _, _, r_sq = _min_image_slab(position[i], position[j], L)
        (r_sq < r_c_sq && !iszero(r_sq)) || continue
        t = Es_Cheb_pair(q[i], q[j], uspara_cheb, sqrt(r_sq))
        energy_short_per_atoms[i] += t / 2
        energy_short_per_atoms[j] += t / 2
    end

    return energy_short_per_atoms
end

"""
    short_energy_Cheb(uspara_cheb, r_c, F0, L, position, q; neighbor_list = nothing)

Short-range (real-space) energy from the Chebyshev interpolant of the
sum-of-Gaussians kernel. `position` is array-of-structs and only indexed
(`p[1]`, `p[2]`, `p[3]`), so any AoS point type works.

`L` is the box `(Lx, Ly, Lz)`; x and y are periodic, z is not. **`r_c` must be
strictly less than `min(Lx, Ly) / 2`** — see `_min_image_slab`.

`neighbor_list` is an optional iterable of candidate pairs, of which only the
first two entries of each element are read; see the note above
`_candidate_pairs`. Omitting it falls back to all i < j pairs.
"""
function short_energy_Cheb(uspara_cheb::ChebPoly{1, T, T}, r_c::T, F0::T, L::NTuple{3, T}, position, q::Vector{T}; neighbor_list = nothing) where{T}

    energy_short = _short_Cheb_sum(_candidate_pairs(_pair_list(neighbor_list), length(q)), uspara_cheb, r_c, L, position, q)

    energy_short += Es_self(q, F0)

    return energy_short / 4π
end

"""
    short_energy_Cheb_per_atom(uspara_cheb, r_c, F0, L, position, q; neighbor_list = nothing)

Per-atom decomposition of [`short_energy_Cheb`](@ref); pair terms are split
evenly between the two partners.
"""
function short_energy_Cheb_per_atom(uspara_cheb::ChebPoly{1, T, T}, r_c::T, F0::T, L::NTuple{3, T}, position, q::Vector{T}; neighbor_list = nothing) where{T}

    energy_short_per_atoms = zeros(T, length(q))

    _short_Cheb_sum_per_atom!(energy_short_per_atoms, _candidate_pairs(_pair_list(neighbor_list), length(q)), uspara_cheb, r_c, L, position, q)

    for i in 1:length(q)
        energy_short_per_atoms[i] += q[i]^2 * F0
    end

    return energy_short_per_atoms ./ 4π
end

function energy_short(interaction::FSSoGInteraction{T}; neighbor_list = nothing) where{T}
    return short_energy_Cheb(interaction.uspara_cheb, interaction.r_c, interaction.F0, interaction.L, interaction.position, interaction.charge; neighbor_list = neighbor_list) / interaction.ϵ
end

function energy_short_per_atom(interaction::FSSoGInteraction{T}; neighbor_list = nothing) where{T}
    return short_energy_Cheb_per_atom(interaction.uspara_cheb, interaction.r_c, interaction.F0, interaction.L, interaction.position, interaction.charge; neighbor_list = neighbor_list) ./ interaction.ϵ
end

function energy_short(interaction::FSSoGThinInteraction{T}; neighbor_list = nothing) where{T}
    return short_energy_Cheb(interaction.uspara_cheb, interaction.r_c, interaction.F0, interaction.L, interaction.position, interaction.charge; neighbor_list = neighbor_list) / interaction.ϵ
end

function energy_short_per_atom(interaction::FSSoGThinInteraction{T}; neighbor_list = nothing) where{T}
    return short_energy_Cheb_per_atom(interaction.uspara_cheb, interaction.r_c, interaction.F0, interaction.L, interaction.position, interaction.charge; neighbor_list = neighbor_list) ./ interaction.ϵ
end

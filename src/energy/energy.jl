# ---------------------------------------------------------------------------
# The framework-free query surface
# ---------------------------------------------------------------------------
#
# Every entry point here takes positions and charges as ARGUMENTS rather than
# reading them from the plan. `FSSoGInteraction` and `FSSoGThinInteraction`
# still carry `position`/`charge` fields, but those are now plan-owned SCRATCH:
# a query scatters the caller's arrays into them and the kernels read the
# scratch. So the caller's data is never mutated, and nothing is allocated per
# call.
#
# `energy` is defined but deliberately NOT exported. Several sibling
# electrostatics packages define a function of the same name, and exporting it
# from all of them would make the bare name ambiguous on
# `using FastSpecSoG, SomeOtherPackage`. Call it as
# `FastSpecSoG.energy(plan, poses, charges)`. The descriptive names
# (`energy_short`, `energy_mid`, `energy_long`, `energy_naive`,
# `energy_per_atom`, ...) are not ambiguous and stay exported.

"""
    _scatter!(interaction, poses, charges)

Copy the caller's positions and charges into the plan's own scratch arrays.
`poses` is array-of-structs and is only indexed (`p[1]`, `p[2]`, `p[3]`), so
`NTuple{3,T}`, `SVector{3,T}` and an MD framework's own point type all work
with no conversion layer at the kernels. Neither `poses` nor `charges` is
modified.
"""
function _scatter!(interaction, poses, charges)
    T = eltype(interaction.charge)
    n_atoms = interaction.n_atoms
    @assert length(poses) == n_atoms "expected $(n_atoms) positions, got $(length(poses))"
    @assert length(charges) == n_atoms "expected $(n_atoms) charges, got $(length(charges))"

    @inbounds for i in 1:n_atoms
        p = poses[i]
        interaction.position[i] = (T(p[1]), T(p[2]), T(p[3]))
        interaction.charge[i] = T(charges[i])
    end

    return interaction
end

"""
    energy_naive(interaction::FSSoG_naive, poses, charges; neighbor_list = nothing)

Total electrostatic energy by direct summation — the internal reference the
fast paths are validated against. `O(N^2 K)`.

`poses` is array-of-structs and is only indexed; nothing is mutated.
`neighbor_list` is an optional iterable of candidate pairs for the short-range
part (only the first two entries of each element are read, and the true
three-dimensional separation is always recomputed from the positions); omitting
it falls back to all i < j pairs.

**`interaction.r_c` must be strictly less than `min(Lx, Ly) / 2`.**
"""
function energy_naive(interaction::FSSoG_naive{T}, poses, charges; neighbor_list = nothing) where{T}

    E_long = long_energy_naive(interaction, poses, charges)
    E_short = short_energy_naive(interaction, poses, charges; neighbor_list = neighbor_list)

    return E_long + E_short
end

"""
    FastSpecSoG.energy(interaction, poses, charges; neighbor_list = nothing)

Total electrostatic energy of the quasi-2D system described by `interaction`,
which may be an [`FSSoGInteraction`](@ref) (cubic geometry: short + middle +
long range), an [`FSSoGThinInteraction`](@ref) (thin slab: short + long range)
or an [`FSSoG_naive`](@ref) (direct summation).

`poses` is array-of-structs and is only indexed (`p[1]`, `p[2]`, `p[3]`);
neither `poses` nor `charges` is mutated — they are scattered into the plan's
own scratch. `neighbor_list` is an optional iterable of candidate pairs for the
short-range part; omitting it falls back to all i < j pairs.

**`interaction.r_c` must be strictly less than `min(Lx, Ly) / 2`.**

Not exported; call it as `FastSpecSoG.energy(...)`.
"""
function energy(interaction::FSSoGInteraction{T}, poses, charges; neighbor_list = nothing) where{T}

    _scatter!(interaction, poses, charges)

    E_long = energy_long(interaction)
    E_mid = energy_mid(interaction)
    E_short = energy_short(interaction; neighbor_list = neighbor_list)

    return E_long + E_mid + E_short
end

function energy(interaction::FSSoGThinInteraction{T}, poses, charges; neighbor_list = nothing) where{T}

    _scatter!(interaction, poses, charges)

    E_long = energy_long(interaction)
    E_short = energy_short(interaction; neighbor_list = neighbor_list)

    return E_long + E_short
end

energy(interaction::FSSoG_naive, poses, charges; neighbor_list = nothing) =
    energy_naive(interaction, poses, charges; neighbor_list = neighbor_list)

"""
    energy_per_atom(interaction, poses, charges; neighbor_list = nothing)

Per-atom decomposition of [`FastSpecSoG.energy`](@ref), in the order of
`poses`. Its sum equals the total energy up to floating-point accumulation
order.
"""
function energy_per_atom(interaction::FSSoGInteraction{T}, poses, charges; neighbor_list = nothing) where{T}

    _scatter!(interaction, poses, charges)

    E_long = energy_long_per_atom(interaction)
    E_mid = energy_mid_per_atom(interaction)
    E_short = energy_short_per_atom(interaction; neighbor_list = neighbor_list)

    return E_long + E_mid + E_short
end

function energy_per_atom(interaction::FSSoGThinInteraction{T}, poses, charges; neighbor_list = nothing) where{T}

    _scatter!(interaction, poses, charges)

    E_long = energy_long_per_atom(interaction)
    E_short = energy_short_per_atom(interaction; neighbor_list = neighbor_list)

    return E_long + E_short
end

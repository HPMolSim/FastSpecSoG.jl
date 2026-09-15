module FastSpecSoGExTinyMDExt

# ExTinyMD adapter for FastSpecSoG. Supplies `ExTinyMD.energy` and nothing else.
#
# ## Why there is no `update_acceleration!`, and no wrapper type
#
# FastSpecSoG is energy-only. It has never computed forces -- there is no
# `force`, no `force!` and no `update_acceleration!` anywhere in the package --
# so it cannot drive an MD run whatever this file does. `simulate!` calls
# `update_acceleration!` on every entry of `sys.interactions` at every step.
#
# That also settles whether the plan types need the dispatcher-wrapper pattern
# that the other packages in this phase use. `MDSys` requires
# `interactions::Vector{T}` with
# `T <: Tuple{ExTinyMD.AbstractInteraction, ExTinyMD.AbstractNeighborFinder}`,
# and a struct's supertype is fixed where the struct is defined: Julia has no
# mechanism for an extension, loaded later and conditionally, to retroactively
# add a supertype to an already-compiled type. `FSSoG_naive`,
# `FSSoGInteraction` and `FSSoGThinInteraction` are defined in this package's
# src/, which does not depend on ExTinyMD at all, so they are not subtypes of
# `ExTinyMD.AbstractInteraction` -- structurally, in every build. A wrapper
# type defined here could carry that supertype, but there would be nothing for
# it to do: with no `update_acceleration!` it would fail on the first step of
# `simulate!`.
#
# Nor is there a working MD integration being given up. Before this change the
# package defined
#
#     ExTinyMD.energy(interaction, neighbor, info::SimulationInfo, atoms::Vector{Atom})
#
# but ExTinyMD's MD loop calls `energy(interaction, neighborfinder, sys, info)`
# (src/MD_core/recorder/energy_logger.jl:51, the only call site) -- a different
# argument order and different types in positions 3 and 4. Those methods could
# therefore never be dispatched by `simulate!` or by `EnergyLogger`; only this
# package's own tests ever called them, directly. They are replaced below by
# methods carrying the signature the MD loop actually uses, which makes this
# adapter uniform with every other one in this phase and with ExTinyMD's own
# src/interactions/electrostatics/adapter.jl.
#
# To be precise about what that does NOT buy: `EnergyLogger` cannot log one
# either. Its constructor requires
# `interactions::Vector{Tuple{T_interaction, T_neighbor}}` with
# `T_interaction <: AbstractInteraction`, and its struct field is typed
# `Vector{Tuple{AbstractInteraction, AbstractNeighborFinder}}`, so the same
# supertype requirement shuts that door too. `ExTinyMD.energy` on a FastSpecSoG
# plan is callable directly and only directly. test/adapter.jl asserts both
# refusals so this comment cannot quietly go stale.
#
# This is the same position ParticleMeshEwald ended in; see
# ExTinyMD.jl/docs/superpowers/specs/2026-09-15-downstream-decoupling-design.md
# sections 4.3a and 5.1.

using FastSpecSoG
using ExTinyMD

"""
    _gather_positions(info) -> Vector{NTuple{3,T}}

Positions in STORAGE-SLOT order. `info.particle_info` is indexed by slot, so
`i` here means "slot i" and matches the charges gathered below.
"""
function _gather_positions(info::ExTinyMD.SimulationInfo{T}) where {T}
    poses = Vector{NTuple{3, T}}(undef, length(info.particle_info))
    @inbounds for i in eachindex(info.particle_info)
        p = info.particle_info[i].position
        poses[i] = (p[1], p[2], p[3])
    end
    return poses
end

"""
    _gather_charges(sys, info) -> Vector{T}

Charges in storage-slot order, honouring ExTinyMD's id/slot indirection:
`sys.atoms` is indexed by particle **id**, `info.particle_info` by storage
**slot**, and the two coincide only in a freshly built `SimulationInfo`.
Reading `sys.atoms[i]` instead of `sys.atoms[info.particle_info[i].id]` is
invisible whenever every particle carries the same charge.
"""
function _gather_charges(sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    charges = Vector{T}(undef, length(info.particle_info))
    @inbounds for i in eachindex(info.particle_info)
        charges[i] = sys.atoms[info.particle_info[i].id].charge
    end
    return charges
end

# A NoNeighborFinder carries no usable list; FastSpecSoG's short-range sums
# fall back to all i < j pairs in that case. Every other finder exposes one,
# and it is treated as candidate pairs only -- the true separation is always
# recomputed from the positions.
_finder_list(::ExTinyMD.NoNeighborFinder) = nothing
_finder_list(f) = f.neighbor_list

const FSSoGPlan = Union{FastSpecSoG.FSSoG_naive,
                        FastSpecSoG.FSSoGInteraction,
                        FastSpecSoG.FSSoGThinInteraction}

"""
    ExTinyMD.energy(interaction, neighborfinder, sys::MDSys, info::SimulationInfo)

Total electrostatic energy of a FastSpecSoG plan, driven from ExTinyMD's MD
types. Refreshes `neighborfinder` if it is due, gathers positions and charges
in storage-slot order (honouring the id/slot indirection), and calls
`FastSpecSoG.energy`.

The plan types cannot be placed in `sys.interactions`, nor in an
`EnergyLogger` -- both require `ExTinyMD.AbstractInteraction`; see the comment
at the top of this file -- so this is callable directly and only directly.
`ExTinyMD.update_acceleration!` is deliberately NOT provided: FastSpecSoG
computes no forces, so there is no way to integrate with it.

Positions and charges are gathered into fresh arrays on each call. That is
deliberate rather than an oversight: this method cannot be on an MD hot path,
since the plans cannot enter `sys.interactions`, so a per-call allocation costs
nothing that matters and avoids a second copy of the plan's own scratch.
"""
function ExTinyMD.energy(interaction::FSSoGPlan, neighborfinder,
                         sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = _gather_positions(info)
    charges = _gather_charges(sys, info)
    return FastSpecSoG.energy(interaction, poses, charges;
                              neighbor_list = _finder_list(neighborfinder))
end

"""
    FastSpecSoG.energy_per_atom(interaction, neighborfinder, sys::MDSys, info::SimulationInfo)

Per-atom decomposition of the energy above, in STORAGE-SLOT order (so entry `i`
belongs to `info.particle_info[i]`, whose particle id is
`info.particle_info[i].id`). Only defined for the two fast plans;
`FSSoG_naive` has no per-atom path.
"""
function FastSpecSoG.energy_per_atom(interaction::Union{FastSpecSoG.FSSoGInteraction,
                                                        FastSpecSoG.FSSoGThinInteraction},
                                     neighborfinder, sys::ExTinyMD.MDSys{T},
                                     info::ExTinyMD.SimulationInfo{T}) where {T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = _gather_positions(info)
    charges = _gather_charges(sys, info)
    return FastSpecSoG.energy_per_atom(interaction, poses, charges;
                                       neighbor_list = _finder_list(neighborfinder))
end

end

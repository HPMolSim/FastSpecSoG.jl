# The ExTinyMD adapter in ext/FastSpecSoGExTinyMDExt.jl.
#
# What this can and cannot test, stated up front so the coverage is not
# overclaimed. FastSpecSoG is energy-only: no `force`, no `force!`, no
# `update_acceleration!`. So there is nothing to drive through `simulate!`, and
# the "adapter under load" test that the other packages in this phase need does
# not apply here -- there is no mass division and no acceleration accumulation
# to get wrong. What IS testable, and is the fault this file exists to catch, is
# the id/slot indirection in the gather.

@testset "adapter: types cannot enter sys.interactions" begin
    # This pins the reasoning in the extension's header comment rather than
    # leaving it as prose. A struct's supertype is fixed where the struct is
    # defined, and these are defined in src/, which has no ExTinyMD.
    @test !(FSSoGInteraction{Float64} <: ExTinyMD.AbstractInteraction)
    @test !(FSSoGThinInteraction{Float64} <: ExTinyMD.AbstractInteraction)
    @test !(FSSoG_naive{Float64} <: ExTinyMD.AbstractInteraction)

    # ... and therefore MDSys rejects them, which is why the extension supplies
    # `ExTinyMD.energy` for direct calls and nothing more.
    L = 20.0
    boundary = Q2dBoundary(L, L, L)
    atoms = [Atom(type = 1, mass = 1.0, charge = 1.0),
             Atom(type = 1, mass = 1.0, charge = -1.0)]
    # r_c = 4.0 < min(Lx, Ly)/2 = 10.0
    plan = FSSoG_naive((L, L, L), 2, 4.0, 2.0, preset = 3)
    @test_throws MethodError MDSys(n_atoms = 2, atoms = atoms, boundary = boundary,
                                   interactions = [(plan, NoNeighborFinder())],
                                   loggers = [TemperatureLogger(10^9; output = false)],
                                   simulator = VerletProcess(dt = 1e-3))

    # `update_acceleration!` is deliberately absent, so an MD run is impossible
    # regardless of the above.
    @test isempty(methods(ExTinyMD.update_acceleration!, (FSSoGInteraction{Float64}, Any, Any, Any)))
    @test isempty(methods(ExTinyMD.update_acceleration!, (FSSoGThinInteraction{Float64}, Any, Any, Any)))
end

# Build an MDSys/SimulationInfo pair with a NoInteraction placeholder (the plans
# cannot go in there, per the testset above).
function _adapter_system(n_atoms, L, Lz, charge_of)
    boundary = Q2dBoundary(L, L, Lz)
    atoms = [Atom(type = 1, mass = 1.0, charge = charge_of(i)) for i in 1:n_atoms]
    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, Lz), boundary;
                          min_r = 1.0, temp = 1.0)
    info.running_step = 1
    sys = MDSys(n_atoms = n_atoms, atoms = atoms, boundary = boundary,
                interactions = [(NoInteraction(), NoNeighborFinder())],
                loggers = [TemperatureLogger(10^9; output = false)],
                simulator = VerletProcess(dt = 1e-3))
    return boundary, atoms, info, sys
end

@testset "adapter: agrees with a direct plan query" begin
    Random.seed!(20260916)
    n_atoms = 40
    L = 50.0
    r_c = 10.0      # r_c = 10.0 < min(Lx, Ly)/2 = 25.0
    boundary, atoms, info, sys = _adapter_system(n_atoms, L, L, i -> isodd(i) ? 2.0 : -2.0)
    finder = CellList3D(info, r_c, boundary, 1)

    # Positions come out of `info` as ExTinyMD `Point{3,Float64}` values and are
    # handed to the core with no conversion: read-only position parameters are
    # untyped and index with p[1]/p[2]/p[3], which is what makes that possible.
    poses = [info.particle_info[i].position for i in 1:n_atoms]
    charges = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]
    @test eltype(poses) <: ExTinyMD.Point

    for plan in (FSSoG_naive((L, L, L), n_atoms, r_c, 4.0, preset = 3),
                 FSSoGInteraction((L, L, L), n_atoms, r_c, 48, 0.5, (128, 128, 128),
                                  (16, 16, 16), 5.0 .* (16, 16, 16), 2, 10, 3,
                                  (32, 32, 32), 48, 32, 32; preset = 3, ϵ = 1.0))
        E_adapter = ExTinyMD.energy(plan, finder, sys, info)
        E_direct = FastSpecSoG.energy(plan, poses, charges;
                                      neighbor_list = finder.neighbor_list)
        # The adapter is a gather plus the same call, so this is exact, not
        # approximate. Any tolerance here would hide a real fault.
        @test E_adapter == E_direct
        @test isfinite(E_adapter)

        # A NoNeighborFinder carries no list; the core then falls back to all
        # i < j pairs, which must give the same energy (different summation
        # order only). Measured relative difference 0.0 for both plans.
        E_nofinder = ExTinyMD.energy(plan, NoNeighborFinder(), sys, info)
        @test E_nofinder ≈ E_adapter rtol = 1e-14
    end
end

@testset "adapter: correct when slot order differs from id order" begin
    # In a freshly built SimulationInfo `particle_info[i].id == i`, so a gather
    # that reads `sys.atoms[i]` instead of `sys.atoms[particle_info[i].id]`
    # looks correct forever. Reverse the slot order so the two differ, and give
    # every id a DISTINCT charge so a mix-up cannot cancel out -- with equal
    # charges the bug is invisible no matter how the slots are permuted.
    Random.seed!(20260917)
    n_atoms = 40
    L = 50.0
    r_c = 10.0      # r_c = 10.0 < min(Lx, Ly)/2 = 25.0

    # Distinct per id, and neutral overall: +-(1 + i/10) in cancelling pairs,
    # so charges run from 1.2 to 5.0 in magnitude. The spread is deliberately
    # wide -- it sets how large an id/slot fault has to be before the assertion
    # below can see it.
    charge_of(i) = (isodd(i) ? 1.0 : -1.0) * (1.0 + (i + isodd(i)) / 10)
    boundary, atoms, info, sys = _adapter_system(n_atoms, L, L, charge_of)
    @test abs(sum(a.charge for a in atoms)) < 1e-12
    @test length(unique(a.charge for a in atoms)) == n_atoms

    reverse!(info.particle_info)
    info.id_dict = Dict(info.particle_info[slot].id => slot for slot in 1:n_atoms)
    # slot order and id order really do differ now
    @test [p.id for p in info.particle_info] == collect(n_atoms:-1:1)
    @test any(i -> info.particle_info[i].id != i, 1:n_atoms)

    finder = CellList3D(info, r_c, boundary, 1)

    plan = FSSoGInteraction((L, L, L), n_atoms, r_c, 48, 0.5, (128, 128, 128),
                            (16, 16, 16), 5.0 .* (16, 16, 16), 2, 10, 3,
                            (32, 32, 32), 48, 32, 32; preset = 3, ϵ = 1.0)

    E_adapter = ExTinyMD.energy(plan, finder, sys, info)

    # The correct gather: positions in slot order, charges by the id stored in
    # that slot.
    poses = [info.particle_info[i].position for i in 1:n_atoms]
    charges_by_id = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]
    E_correct = FastSpecSoG.energy(plan, poses, charges_by_id;
                                   neighbor_list = finder.neighbor_list)
    @test E_adapter == E_correct

    # The fault this testset exists for: charges read by SLOT instead of by id.
    # Asserting that the adapter differs from this is what makes the test above
    # non-vacuous -- without it, both would pass an adapter that ignored
    # particle_info[i].id entirely.
    charges_by_slot = [atoms[i].charge for i in 1:n_atoms]
    @test charges_by_slot != charges_by_id
    E_wrong = FastSpecSoG.energy(plan, poses, charges_by_slot;
                                 neighbor_list = finder.neighbor_list)
    # Measured |E_correct - E_wrong| = 6.066e-2 on this configuration, against
    # a total of about -0.56, i.e. an 11% error. The bound of 1e-3 is 61x below
    # the measured value; it is a "this must not be small" assertion, not an
    # accuracy one. (With the charges spread only over 1.02-1.40 instead of
    # 1.2-5.0 the same fault measured 2.368e-3, which is why the spread is
    # wide: a narrow one makes an id/slot fault nearly invisible.)
    @test abs(E_correct - E_wrong) > 1e-3

    # The per-atom adapter is in slot order too, and sums to the total.
    # Measured residual 4.219e-15 against a 1e-13 bound (24x).
    E_pa = FastSpecSoG.energy_per_atom(plan, finder, sys, info)
    @test length(E_pa) == n_atoms
    @test all(isfinite, E_pa)
    @test abs(sum(E_pa) - E_adapter) < 1e-13
end

@testset "adapter: thin plan" begin
    Random.seed!(20260918)
    n_atoms = 40
    L = 100.0
    Lz = 1.0
    r_c = 10.0      # r_c = 10.0 < min(Lx, Ly)/2 = 50.0
    boundary, atoms, info, sys = _adapter_system(n_atoms, L, Lz, i -> isodd(i) ? 2.0 : -2.0)
    finder = CellList3D(info, r_c, boundary, 1)

    plan = FSSoGThinInteraction((L, L, Lz), n_atoms, r_c, 48, 0.5, (128, 128), 32,
                                (16, 16), 5.0 .* (16, 16), 16, 24, 32, 32;
                                preset = 3, ϵ = 1.0)

    poses = [info.particle_info[i].position for i in 1:n_atoms]
    charges = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]

    E_adapter = ExTinyMD.energy(plan, finder, sys, info)
    @test isfinite(E_adapter)
    @test E_adapter == FastSpecSoG.energy(plan, poses, charges;
                                          neighbor_list = finder.neighbor_list)

    # Measured per-atom sum residual 0.0; bound 1e-13.
    E_pa = FastSpecSoG.energy_per_atom(plan, finder, sys, info)
    @test abs(sum(E_pa) - E_adapter) < 1e-13
end

@testset "adapter: neither simulate! nor EnergyLogger can reach it" begin
    # An honest statement of the adapter's reach, asserted rather than assumed.
    #
    # The extension gives `ExTinyMD.energy` the signature ExTinyMD 0.3's MD loop
    # uses, `(interaction, neighborfinder, sys::MDSys, info::SimulationInfo)`,
    # replacing the old `(interaction, neighbor, info, atoms)` which had an
    # argument order no ExTinyMD code path ever calls. That makes it uniform
    # with every other adapter in this phase and with ExTinyMD's own
    # src/interactions/electrostatics/adapter.jl.
    #
    # It does NOT make the plans usable from `simulate!` or from
    # `EnergyLogger`: both require `AbstractInteraction`, which a type defined
    # in src/ (where ExTinyMD does not exist) can never be. `MDSys` is checked
    # above; `EnergyLogger` is checked here, because its own docstring
    # ("Records ... the energy of each of `interactions`") makes it the
    # plausible place to assume otherwise. Its constructor requires
    # `interactions::Vector{Tuple{T_interaction, T_neighbor}}` with
    # `T_interaction <: AbstractInteraction`, and its struct field is typed
    # `Vector{Tuple{AbstractInteraction, AbstractNeighborFinder}}`, so there is
    # no way in either.
    #
    # Conclusion: `ExTinyMD.energy(plan, finder, sys, info)` is callable
    # directly and only directly. That is ParticleMeshEwald's position (design
    # doc sections 4.3a and 5.1) and is fine here, because FastSpecSoG computes
    # no forces and so could not drive an MD run in any case.
    Random.seed!(20260919)
    n_atoms = 20
    L = 50.0
    r_c = 10.0      # r_c = 10.0 < min(Lx, Ly)/2 = 25.0
    boundary, atoms, info, sys = _adapter_system(n_atoms, L, L, i -> isodd(i) ? 2.0 : -2.0)
    finder = CellList3D(info, r_c, boundary, 1)
    plan = FSSoG_naive((L, L, L), n_atoms, r_c, 4.0, preset = 3)

    @test_throws MethodError EnergyLogger(1; interactions = [(plan, finder)],
                                          energy_names = ["fssog"], output = false)

    # ... while the direct call works and is what the README documents.
    @test isfinite(ExTinyMD.energy(plan, finder, sys, info))
end

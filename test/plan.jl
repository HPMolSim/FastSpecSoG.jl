# The framework-free query API.
#
# Nothing in this file constructs, names or touches an ExTinyMD type: the plans
# are built from a box and an atom count, and queried with plain arrays. That
# is the contract this phase exists to establish. It is *not* proof that
# ExTinyMD is unnecessary -- test/runtests.jl loads ExTinyMD for the adapter
# tests, so it is already in this process. test/standalone.jl proves that part
# in a fresh subprocess.

@testset "plan API: slab minimum image" begin
    L = (50.0, 50.0, 50.0)

    # The property the whole short-range sum rests on: the returned squared
    # distance is THREE-DIMENSIONAL. An in-plane-only helper -- the substitution
    # the phase plan warns against, and the one that made ExTinyMD's own
    # Ewald2D return +0.0238 against a true -0.1539 before Phase 1 -- would
    # return 0.0 for this pair, because the two particles share an (x, y)
    # column.
    _, _, r_sq = FastSpecSoG._min_image_slab((1.0, 2.0, 0.0), (1.0, 2.0, 3.0), L)
    @test r_sq == 9.0

    # x and y wrap to the nearest image ...
    _, _, r_wrap = FastSpecSoG._min_image_slab((1.0, 1.0, 1.0), (49.0, 1.0, 1.0), L)
    @test r_wrap == 4.0          # an unwrapped difference would give 48^2 = 2304
    _, _, r_wrap_y = FastSpecSoG._min_image_slab((1.0, 1.0, 1.0), (1.0, 49.0, 1.0), L)
    @test r_wrap_y == 4.0

    # ... and z, the slab axis, does not.
    _, _, r_z = FastSpecSoG._min_image_slab((1.0, 1.0, 1.0), (1.0, 1.0, 49.0), L)
    @test r_z == 2304.0          # a wrapping z would give 4.0

    # Returned coordinates are NTuple, not SVector: this package has no
    # StaticArrays dependency and the coordinates are only ever indexed.
    ci, cj, _ = FastSpecSoG._min_image_slab((1.0, 1.0, 1.0), (49.0, 1.0, 1.0), L)
    @test ci isa NTuple{3, Float64}
    @test cj isa NTuple{3, Float64}
    # coord_i is pos_i's nearest in-plane image of pos_j, coord_j is pos_j itself
    @test ci == (51.0, 1.0, 1.0)
    @test cj == (49.0, 1.0, 1.0)

    # Inputs need only support indexing.
    _, _, r_idx = FastSpecSoG._min_image_slab(IdxOnlyPoint((1.0, 1.0, 1.0)),
                                              IdxOnlyPoint((49.0, 1.0, 1.0)), L)
    @test r_idx == 4.0
end

@testset "plan API: short range sees the out-of-plane separation" begin
    # An energy-level version of the test above, which is the one that would
    # catch an in-plane substitution buried inside the sum rather than in the
    # helper. Two charges share an (x, y) column and are 2.0 apart in z, well
    # inside the cutoff. An in-plane distance would make r_sq zero, the pair
    # would be rejected as coincident, and only the self term would survive.
    L = (10.0, 10.0, 10.0)
    r_c = 4.0      # r_c = 4.0 < min(Lx, Ly)/2 = 5.0
    poses = [(3.0, 3.0, 2.0), (3.0, 3.0, 4.0)]
    charges = [1.0, -1.0]

    plan = FSSoG_naive(L, 2, r_c, 3.0, preset = 3)
    E = short_energy_naive(plan, poses, charges)
    E_self = FastSpecSoG.Es_naive_self(charges, plan) / (4π * plan.ϵ)

    # The pair term is q_i q_j (1/r - U(r)) / 4pi at r = 2.0, negative here.
    E_pair = FastSpecSoG.Es_naive_pair(1.0, -1.0, plan.uspara, 4.0) / (4π * plan.ϵ)
    @test E_pair < 0
    @test E ≈ E_self + E_pair rtol = 1e-14
    # and it is not merely the self term
    @test !isapprox(E, E_self; rtol = 1e-6)
end

@testset "plan API: cube geometry against the naive reference" begin
    Random.seed!(20260915)
    n_atoms = 100
    L = (50.0, 50.0, 50.0)
    r_c = 10.0     # r_c = 10.0 < min(Lx, Ly)/2 = 25.0
    poses = _random_poses(n_atoms, L; min_r = 1.0)
    charges = [isodd(i) ? 2.0 : -2.0 for i in 1:n_atoms]

    poses_ref = deepcopy(poses)
    charges_ref = deepcopy(charges)

    naive = FSSoG_naive(L, n_atoms, r_c, 4.0, preset = 3)
    E_naive = energy_naive(naive, poses, charges)

    plan = FSSoGInteraction(L, n_atoms, r_c, 48, 0.5, (128, 128, 128), (16, 16, 16),
                            5.0 .* (16, 16, 16), 2, 10, 3, (32, 32, 32), 48, 32, 32;
                            preset = 3, ϵ = 1.0)
    E = FastSpecSoG.energy(plan, poses, charges)

    # Measured |E - E_naive| for this configuration: 1.2641e-10. Bound 1e-9,
    # a margin of 7.9x. This is the fast method against its own direct sum, so
    # the residual is the method's truncation error, not round-off; the bound is
    # set just above it so a regression in any of the three ranges breaks it.
    # (For comparison, test/energy.jl holds the same comparison to 1e-6.)
    @test abs(E - E_naive) < 1e-9

    # The long+middle ranges together are the reciprocal-space sum that
    # long_energy_naive computes directly. Measured difference 1.2641e-10;
    # same 1e-9 bound (7.9x) and the same reasoning.
    E_long_naive = long_energy_naive(naive, poses, charges)
    @test abs((energy_long(plan) + energy_mid(plan)) - E_long_naive) < 1e-9

    # The query must not touch the caller's arrays.
    @test poses == poses_ref
    @test charges == charges_ref

    # Per-atom decomposition sums to the total. Measured residual 5.00e-16,
    # bound 1e-13 (200x). The margin is deliberately wide: the failure this
    # guards against is a wrong decomposition -- a dropped or double-counted
    # range, or half a pair term assigned to one partner -- which is an O(1)
    # error, not a round-off drift, so a tight bound would buy nothing and
    # would break on a different thread count or machine.
    E_per_atom = energy_per_atom(plan, poses, charges)
    @test length(E_per_atom) == n_atoms
    @test all(isfinite, E_per_atom)
    @test abs(sum(E_per_atom) - E) < 1e-13
    @test poses == poses_ref
    @test charges == charges_ref

    # Supplying a candidate-pair list must give the same answer as the
    # all-pairs fallback: the list is filtered on a recomputed distance, so
    # only the summation order changes.
    nb = [(i, j, 0.0) for i in 1:n_atoms - 1 for j in i + 1:n_atoms]
    @test FastSpecSoG.energy(plan, poses, charges; neighbor_list = nb) ≈ E rtol = 1e-14

    # Positions need only support indexing.
    poses_idx = [IdxOnlyPoint(p) for p in poses]
    @test FastSpecSoG.energy(plan, poses_idx, charges) ≈ E rtol = 1e-14
    @test energy_naive(naive, poses_idx, charges) ≈ E_naive rtol = 1e-14

    # Length mismatches are caught rather than read out of bounds.
    @test_throws AssertionError FastSpecSoG.energy(plan, poses[1:end - 1], charges)
    @test_throws AssertionError FastSpecSoG.energy(plan, poses, charges[1:end - 1])
end

@testset "plan API: thin geometry against the naive reference" begin
    Random.seed!(20260915)
    n_atoms = 100
    L = (100.0, 100.0, 1.0)
    r_c = 10.0     # r_c = 10.0 < min(Lx, Ly)/2 = 50.0
    poses = _random_poses(n_atoms, L; min_r = 1.0)
    charges = [isodd(i) ? 2.0 : -2.0 for i in 1:n_atoms]

    poses_ref = deepcopy(poses)
    charges_ref = deepcopy(charges)

    naive = FSSoG_naive(L, n_atoms, r_c, 3.0, preset = 3)
    E_naive = energy_naive(naive, poses, charges)

    plan = FSSoGThinInteraction(L, n_atoms, r_c, 48, 0.5, (128, 128), 32, (16, 16),
                                5.0 .* (16, 16), 16, 24, 32, 32; preset = 3, ϵ = 1.0)
    E = FastSpecSoG.energy(plan, poses, charges)

    # Measured |E - E_naive| for this configuration: 2.8223e-10. Bound 3e-9,
    # a margin of 10.6x. (test/energy.jl holds the same comparison to 1e-6.)
    @test abs(E - E_naive) < 3e-9

    # Measured per-atom sum residual 1.554e-15, bound 1e-13 (64x); see the
    # cube testset for why this margin is wide on purpose.
    E_per_atom = energy_per_atom(plan, poses, charges)
    @test all(isfinite, E_per_atom)
    @test abs(sum(E_per_atom) - E) < 1e-13

    @test poses == poses_ref
    @test charges == charges_ref

    poses_idx = [IdxOnlyPoint(p) for p in poses]
    @test FastSpecSoG.energy(plan, poses_idx, charges) ≈ E rtol = 1e-14
end

@testset "plan API: degenerate separations stay finite" begin
    # Replacing `position_check3D`'s all-zero sentinel with an explicit cutoff
    # test keeps the pairs the sentinel was incidentally discarding. Both
    # `Es_Cheb_pair` and `Es_naive_pair` divide by r, so any r == 0 pair that
    # survives the guard produces Inf -- and NaN once a zero charge multiplies
    # it. In QuasiEwald that exact conversion shipped a 0/0 that fired for any
    # lattice initialisation and was caught only in final review.
    L = (50.0, 50.0, 50.0)
    r_c = 10.0     # r_c = 10.0 < min(Lx, Ly)/2 = 25.0
    n_atoms = 20

    Random.seed!(777)
    # Particles 1 and 2 sit at exactly the same site; 3 and 4 are separated by
    # exactly Lx in x with identical y and z, so their minimum image is also
    # r == 0. The remaining 16 are placed at least min_r from all of them, so
    # every pair that is NOT degenerate stays inside the Chebyshev
    # interpolant's domain and a domain error cannot mask a NaN.
    degenerate = [(10.0, 10.0, 10.0), (10.0, 10.0, 10.0),
                  (0.0, 12.5, 7.5), (L[1], 12.5, 7.5)]
    poses = _random_poses(n_atoms - 4, L; min_r = 1.0, initial = degenerate)
    @test length(poses) == n_atoms
    @test _slab_r2(poses[1], poses[2], L) == 0.0
    @test _slab_r2(poses[3], poses[4], L) == 0.0
    # The coincident pair carries zero charge, so an unguarded 1/r gives
    # 0 * Inf = NaN rather than an obvious Inf; the exact-Lx pair carries
    # non-zero charge, so it gives Inf. Both failure modes are covered, and the
    # charges stay neutral overall.
    charges = vcat([0.0, 0.0, 2.0, -2.0],
                   [isodd(i) ? 2.0 : -2.0 for i in 5:n_atoms])
    @test sum(charges) == 0.0

    naive = FSSoG_naive(L, n_atoms, r_c, 4.0, preset = 3)
    E_naive = energy_naive(naive, poses, charges)
    @test isfinite(E_naive)

    uspara_cheb, F0 = Es_Cheb_precompute(3, 0.5, r_c, 64)
    E_cheb = short_energy_Cheb(uspara_cheb, r_c, F0, L, poses, charges)
    @test isfinite(E_cheb)
    E_cheb_pa = FastSpecSoG.short_energy_Cheb_per_atom(uspara_cheb, r_c, F0, L, poses, charges)
    @test all(isfinite, E_cheb_pa)

    plan = FSSoGInteraction(L, n_atoms, r_c, 48, 0.5, (128, 128, 128), (16, 16, 16),
                            5.0 .* (16, 16, 16), 2, 10, 3, (32, 32, 32), 48, 32, 32;
                            preset = 3, ϵ = 1.0)
    E = FastSpecSoG.energy(plan, poses, charges)
    @test isfinite(E)
    @test all(isfinite, energy_per_atom(plan, poses, charges))

    # Skipping the r == 0 pairs is also the pre-existing behaviour, so the two
    # paths still agree. Measured |E - E_naive| = 6.288e-11; bound 1e-9, a
    # margin of 15.9x.
    @test abs(E - E_naive) < 1e-9
end

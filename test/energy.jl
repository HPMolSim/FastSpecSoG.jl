# FastSpecSoG against an independent Ewald2D reference.
#
# The reference used to be EwaldSummations.jl, which is unmaintained and pinned
# to CellListMap 0.9 and therefore cannot co-resolve with ExTinyMD 0.3. It is
# replaced by ExTinyMD's own `Ewald2D`, the Phase 1 implementation validated to
# ~5e-9 α-independence and against the image-charge method to 4.7e-8.
#
# The API change is not only a rename. The old call was POSITIONAL and took s
# before α:
#
#     Ewald2DInteraction(n_atoms, 5.0, 0.25, (L, L, L), ϵ = 1.0)   # (n, s, α, L)
#
# The new one takes both as KEYWORDS:
#
#     Ewald2D(n_atoms, (L, L, L); α = 0.25, s = 5.0, ϵ = 1.0)
#
# Getting them the wrong way round does not error, and the total Ewald energy is
# independent of α, so the total would still be right while the splitting was
# wrong. Every test below therefore pins the SPLIT -- the interaction's own
# α/r_c/k_c fields, and the short and long parts separately against
# `ewald2d_per_atom`, which recomputes r_c and k_c from literal α and s rather
# than reading them off the interaction.

@testset "energy cube" begin

    @info "testing the energy cube"

    n_atoms = 100
    L = 50.0
    Lt = (L, L, L)

    α = 0.25
    s = 5.0
    # Ewald r_c = s/α = 20.0 < min(Lx, Ly)/2 = 25.0
    ewald = Ewald2D(n_atoms, Lt; α = α, s = s, ϵ = 1.0)

    # The α/s swap guard. Verified by construction: with the pair transposed
    # (α = 5.0, s = 0.25) the interaction reports α = 5.0 and r_c = 0.05, so
    # the α and r_c assertions fail. NOTE that k_c does NOT catch a swap --
    # k_c = 2αs is symmetric in the two, and measures 2.5 either way -- so the
    # k_c assertions below are a cutoff-formula check, not part of the guard.
    @test ewald.short.α == α
    @test ewald.long.α == α
    @test ewald.short.r_c ≈ s / α
    @test ewald.short.k_c ≈ 2 * α * s
    @test ewald.long.k_c ≈ 2 * α * s

    Random.seed!(20250915)
    poses = _random_poses(n_atoms, Lt; min_r = 1.0)
    charges = [isodd(i) ? 2.0 : -2.0 for i in 1:n_atoms]
    @test sum(charges) == 0.0

    energy_ewald = coulomb_energy(ewald, poses, charges)
    E_s = short_energy(ewald.short, poses, charges)
    E_l = long_energy(ewald.long, poses, charges)
    @test E_s + E_l ≈ energy_ewald rtol = 1e-14

    # The independent per-atom reference, built from literal α and s. Its short
    # and long parts are asserted SEPARATELY, which is what a swapped α/s pair
    # would break; the total alone would not catch it.
    ewald_s_pa, ewald_l_pa = ewald2d_per_atom(poses, charges, Lt; α = α, s = s, ϵ = 1.0)
    energy_ewald_per_atom = ewald_s_pa .+ ewald_l_pa
    # Measured: |sum(ewald_s_pa) - E_s| = 3.553e-15,
    # |sum(ewald_l_pa) - E_l| = 8.88e-16, total 4.00e-15. These are the same
    # sums in a different order, so they are round-off bounds; 1e-12 gives
    # margins of 281x, 1126x and 250x. They are also the assertions that pin
    # the split: the reference's α and s come from the literals above, so a
    # transposed pair at the `Ewald2D` call would make the short and long parts
    # disagree by O(1) (measured below: 4.58 at twice this α) while the total
    # stayed correct.
    @test abs(sum(ewald_s_pa) - E_s) < 1e-12
    @test abs(sum(ewald_l_pa) - E_l) < 1e-12
    @test abs(sum(energy_ewald_per_atom) - energy_ewald) < 1e-12

    # The total is α-independent; the split is not. Doubling α must leave the
    # total alone and move both parts substantially, which is what makes the
    # split assertions above meaningful rather than tautological.
    # Ewald r_c = s/α = 10.0 < 25.0 at this α as well.
    ewald2 = Ewald2D(n_atoms, Lt; α = 2α, s = s, ϵ = 1.0)
    energy_ewald2 = coulomb_energy(ewald2, poses, charges)
    # Measured |E(α) - E(2α)| = 1.298e-12 on this configuration; bound 1e-10,
    # a margin of 77x. Measured movement of each part: 4.578 for both the short
    # and the long, against the lower bound of 1e-3 below.
    @test abs(energy_ewald - energy_ewald2) < 1e-10
    @test abs(short_energy(ewald2.short, poses, charges) - E_s) > 1e-3
    @test abs(long_energy(ewald2.long, poses, charges) - E_l) > 1e-3

    for r_c in [10.0, 15.0]
        # r_c = 10.0 and 15.0, both < min(Lx, Ly)/2 = 25.0
        fssog_naive = FSSoG_naive(Lt, n_atoms, r_c, 4.0, preset = 3)
        energy_sog_naive = energy_naive(fssog_naive, poses, charges)

        N_real = (128, 128, 128)
        w = (16, 16, 16)
        β = 5.0 .* w
        extra_pad_ratio = 2
        cheb_order = 10
        preset = 3
        M_mid = 3

        N_grid = (32, 32, 32)
        Q = 48
        R_z0 = 32
        Q_0 = 32

        fssog_interaction = FSSoGInteraction(Lt, n_atoms, r_c, Q, 0.5, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Q, R_z0, Q_0; preset = preset, ϵ = 1.0)

        energy_sog = FastSpecSoG.energy(fssog_interaction, poses, charges)

        # Bounds unchanged from before the EwaldSummations swap. Measured
        # against the new Ewald2D reference:
        #   r_c = 10: 5.903e-5, 5.903e-5, 5.471e-10
        #   r_c = 15: 6.095e-5, 6.121e-5, 2.645e-7
        # so the margins on the 1e-4 bounds are only ~1.6-1.7x. That is
        # pre-existing: the discrepancy is the sum-of-Gaussians truncation at
        # these grid parameters, and the FastSpecSoG side of the comparison is
        # bit-identical to its pre-Phase-3d values. The bound was NOT loosened
        # to accommodate the new reference. The 1e-6 bound has margins of
        # 1828x and 3.8x.
        @test abs(energy_ewald - energy_sog_naive) < 1e-4
        @test abs(energy_ewald - energy_sog) < 1e-4
        @test abs(energy_sog_naive - energy_sog) < 1e-6

        energy_sog_per_atom = energy_per_atom(fssog_interaction, poses, charges)

        # Bound unchanged. Measured maximum over i in 1:10 -- r_c = 10:
        # 1.582e-6; r_c = 15: 1.904e-6. Margins 63x and 53x.
        for i in 1:10
            @test abs(energy_sog_per_atom[i] - energy_ewald_per_atom[i]) < 1e-4
        end
    end
end

@testset "energy thin" begin

    @info "testing the energy thin"

    n_atoms = 100
    Lx = 100.0
    Ly = 100.0
    Lz = 1.0
    Lt = (Lx, Ly, Lz)

    α = 0.2
    s = 5.0
    # Ewald r_c = s/α = 25.0 < min(Lx, Ly)/2 = 50.0
    ewald = Ewald2D(n_atoms, Lt; α = α, s = s, ϵ = 1.0)

    # Swap guard, as in the cube testset: α and r_c catch a transposed pair,
    # k_c does not.
    @test ewald.short.α == α
    @test ewald.long.α == α
    @test ewald.short.r_c ≈ s / α
    @test ewald.long.k_c ≈ 2 * α * s

    Random.seed!(20250915)
    poses = _random_poses(n_atoms, Lt; min_r = 1.0)
    charges = [isodd(i) ? 2.0 : -2.0 for i in 1:n_atoms]
    @test sum(charges) == 0.0

    energy_ewald = coulomb_energy(ewald, poses, charges)
    E_s = short_energy(ewald.short, poses, charges)
    E_l = long_energy(ewald.long, poses, charges)
    @test E_s + E_l ≈ energy_ewald rtol = 1e-14

    ewald_s_pa, ewald_l_pa = ewald2d_per_atom(poses, charges, Lt; α = α, s = s, ϵ = 1.0)
    energy_ewald_per_atom = ewald_s_pa .+ ewald_l_pa
    # Measured: |sum(ewald_s_pa) - E_s| = 7.550e-15,
    # |sum(ewald_l_pa) - E_l| = 4.885e-15, total 2.442e-15; bound 1e-12 as
    # above (margins 132x, 205x, 409x).
    @test abs(sum(ewald_s_pa) - E_s) < 1e-12
    @test abs(sum(ewald_l_pa) - E_l) < 1e-12
    @test abs(sum(energy_ewald_per_atom) - energy_ewald) < 1e-12

    # Ewald r_c = s/α = 12.5 < 50.0 at 2α as well.
    ewald2 = Ewald2D(n_atoms, Lt; α = 2α, s = s, ϵ = 1.0)
    # Measured |E(α) - E(2α)| = 9.437e-13; bound 1e-10, a margin of 106x.
    # Measured movement of each part: 3.614.
    @test abs(energy_ewald - coulomb_energy(ewald2, poses, charges)) < 1e-10
    @test abs(short_energy(ewald2.short, poses, charges) - E_s) > 1e-3
    @test abs(long_energy(ewald2.long, poses, charges) - E_l) > 1e-3

    for r_c in [10.0, 15.0]
        # r_c = 10.0 and 15.0, both < min(Lx, Ly)/2 = 50.0
        N_real = (128, 128)
        R_z = 32
        w = (16, 16)
        β = 5.0 .* w
        cheb_order = 16
        preset = 3
        Q = 48
        Q_0 = 32
        R_z0 = 32
        Taylor_Q = 24

        fssog_interaction = FSSoGThinInteraction(Lt, n_atoms, r_c, Q, 0.5, N_real, R_z, w, β, cheb_order, Taylor_Q, R_z0, Q_0; preset = preset, ϵ = 1.0)

        energy_sog = FastSpecSoG.energy(fssog_interaction, poses, charges)

        fssog_naive = FSSoG_naive(Lt, n_atoms, r_c, 3.0, preset = 3)
        energy_sog_naive = energy_naive(fssog_naive, poses, charges)

        # Bounds unchanged. Measured: r_c = 10: 2.057e-5, 2.057e-5, 4.463e-11;
        # r_c = 15: 1.730e-5, 1.730e-5, 1.757e-10. Margins ~49-58x on the 1e-3
        # bounds.
        @test abs(energy_ewald - energy_sog) < 1e-3
        @test abs(energy_ewald - energy_sog_naive) < 1e-3
        @test abs(energy_sog - energy_sog_naive) < 1e-6

        energy_sog_per_atom = energy_per_atom(fssog_interaction, poses, charges)
            for i in 1:10
            @test abs(energy_sog_per_atom[i] - energy_ewald_per_atom[i]) < 1e-4
        end
    end
end

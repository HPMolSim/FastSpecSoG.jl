@testset "energy cube" begin

    @info "testing the energy cube"

    n_atoms = 100
    L = 50.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 2.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 2.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 5.0, 0.25, (L, L, L), ϵ = 1.0)
    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)

    p = [info.particle_info[i].position for i in 1:n_atoms]
    q = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]
    energy_ewald_s = Ewald2D_short_energy_N(10, Ewald2D_interaction, p, q)
    energy_ewald_l = Ewald2D_long_energy_N(10, Ewald2D_interaction, p, q)
    energy_ewald_per_atom = energy_ewald_s + energy_ewald_l

    for r_c in [10.0, 15.0]
        fssog_naive = FSSoG_naive((L, L, L), n_atoms, r_c, 4.0, preset = 3)
        fssog_neighbor = CellList3D(info, fssog_naive.r_c, boundary, 1)
        energy_sog_naive = energy_naive(fssog_naive, fssog_neighbor, info, atoms)

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

        fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Q, 0.5, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Q, R_z0, Q_0; preset = preset, ϵ = 1.0)

        fssog_neighbor = CellList3D(info, fssog_interaction.r_c, fssog_interaction.boundary, 1)
        energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

        @test abs(energy_ewald - energy_sog_naive) < 1e-4
        @test abs(energy_ewald - energy_sog) < 1e-4
        @test abs(energy_sog_naive - energy_sog) < 1e-6
        
        energy_sog_per_atom = energy_per_atom(fssog_interaction, fssog_neighbor, info, atoms)

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
    boundary = ExTinyMD.Q2dBoundary(Lx, Ly, Lz)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 2.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 2.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, Lx, 0.0, Ly, 0.0, Lz), boundary; min_r = 1.0, temp = 1.0)

    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 5.0, 0.2, (Lx, Ly, Lz), ϵ = 1.0)
    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)

    p = [info.particle_info[i].position for i in 1:n_atoms]
    q = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]
    energy_ewald_s = Ewald2D_short_energy_N(10, Ewald2D_interaction, p, q)
    energy_ewald_l = Ewald2D_long_energy_N(10, Ewald2D_interaction, p, q)
    energy_ewald_per_atom = energy_ewald_s + energy_ewald_l

    for r_c in [10.0, 15.0]
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

        fssog_interaction = FSSoGThinInteraction((Lx, Ly, Lz), n_atoms, r_c, Q, 0.5, N_real, R_z, w, β, cheb_order, Taylor_Q, R_z0, Q_0; preset = preset, ϵ = 1.0)

        fssog_neighbor = CellList3D(info, fssog_interaction.r_c, fssog_interaction.boundary, 1)
        energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)

        fssog_naive = FSSoG_naive((Lx, Ly, Lz), n_atoms, r_c, 3.0, preset = 3)
        fssog_neighbor = CellList3D(info, fssog_naive.r_c, boundary, 1)
        energy_sog_naive = energy_naive(fssog_naive, fssog_neighbor, info, atoms)

        @test abs(energy_ewald - energy_sog) < 1e-3
        @test abs(energy_ewald - energy_sog_naive) < 1e-3
        @test abs(energy_sog - energy_sog_naive) < 1e-6

        energy_sog_per_atom = energy_per_atom(fssog_interaction, fssog_neighbor, info, atoms)
        for i in 1:10
            @test abs(energy_sog_per_atom[i] - energy_ewald_per_atom[i]) < 1e-4
        end
    end
end
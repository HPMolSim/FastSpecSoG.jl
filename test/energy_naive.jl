@testset "test energy naive" begin
    n_atoms = 100
    L = 100.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    
    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 3.0, 0.2, (L, L, L), ϵ = 1.0)
    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)

    # coords = [p_info.position for p_info in info.particle_info]
    # charge = [atoms[p_info.id].charge for p_info in info.particle_info]
    FastSpecSoG_interaction = FastSpecSOGInteraction((L, L, L), n_atoms)
    FastSpecSoG_neighbor = CellList3D(info, FastSpecSoG_interaction.r_c, boundary, 1)
    energy_sog = energy_naive(FastSpecSoG_interaction, FastSpecSoG_neighbor, info, atoms)

    @show energy_ewald, energy_sog, energy_ewald - energy_sog
    @test abs(energy_ewald - energy_sog) < 1e-4
end
@testset "test energy naive" begin

    @info "testing the energy naive"

    n_atoms = 100
    L = 100.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 2.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 2.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    
    Ewald2D_interaction = Ewald2DInteraction(n_atoms, 5.0, 0.2, (L, L, L), ϵ = 1.0)
    Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
    energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)

    FastSpecSoG_interaction = FSSoG_naive((L, L, L), n_atoms, 10.0, 3.0, preset = 5)
    FastSpecSoG_neighbor = CellList3D(info, FastSpecSoG_interaction.r_c, boundary, 1)
    energy_sog = energy_naive(FastSpecSoG_interaction, FastSpecSoG_neighbor, info, atoms)

    @test abs(energy_ewald - energy_sog) < 1e-4
end
@testset "Cheb approx of USeries_direct" begin
    @info "testing the Cheb approx of USeries_direct"
    for preset in 1:5
        @testset "preset = $preset" begin
            uspara = USeriesPara(preset)
            r_min = 0.5
            r_max = 10.0
            Q = 64
            uspara_cheb, F0 = Es_Cheb_precompute(preset, r_min, r_max, Q)
            x = [0.5:0.01:10.0...]
            f_direct = r -> U_series(r, uspara)
            f_cheb = r -> uspara_cheb(r)
            @test maximum(abs.(f_direct.(x) .- f_cheb.(x))) < 1e-12
        end
    end
end

@testset "Es_Cheb vs Es_direct" begin
    @info "testing the Cheb approx of the short range energy"
    n_atoms = 1000
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
    position = [info.particle_info[i].position.coo for i in 1:n_atoms]
    charge = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]

    for preset in 1:5
        @testset "preset = $preset" begin
            r_min = 0.5
            r_max = 10.0
            Q = 64
            uspara_cheb, F0 = Es_Cheb_precompute(preset, r_min, r_max, Q)
            interaction = FSSoG_naive((L, L, L), n_atoms, r_max, 3.0, preset = preset)
            neighbor = CellList3D(info, interaction.r_c, boundary, 1)

            Es_naive = short_energy_naive(interaction, neighbor, position, charge)
            Es_Cheb = short_energy_Cheb(uspara_cheb, interaction.r_c, F0, boundary, neighbor, position, charge)
            @test isapprox(Es_naive, Es_Cheb; atol = 1e-12)
        end
    end
end
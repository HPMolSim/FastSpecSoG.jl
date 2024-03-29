@testset "energy_long long range" begin
    @info "testing the long range part of the energy"
    n_atoms = 32
    L = (100.0, 100.0, 100.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    N_grid = (32, 32, 32)
    uspara = USeriesPara(2)
    soepara = SoePara16()
    M_mid = 5

    k_x, k_y, r_z, H_r, H_c, phase_x, phase_y, sort_z, z = long_paras_gen(L, N_grid, n_atoms)
    cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, 64)

    E_long_loop_k = energy_long_loop_k(qs, poses, L, M_mid, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, uspara)
    E_long_cheb_k = energy_long_cheb_k(qs, poses, L, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, cheb_mat)

    E_long_0 = energy_long_0(qs, poses, L, M_mid, z, sort_z, uspara, soepara)

    @info "running the direct summation for the long range part of the energy"
    # using the direct summation
    E_direct_k = long_energy_us_k(qs, poses, 50, L, uspara, M_mid + 1, length(uspara.sw))
    E_direct_0 = long_energy_us_0(qs, poses, L, uspara, M_mid + 1, length(uspara.sw))

    @test isapprox(E_long_0, E_direct_0, atol = 1e-4)
    @test isapprox(E_long_loop_k, E_long_cheb_k)
    @test isapprox(E_long_loop_k, E_direct_k, atol=1e-8)
    @test isapprox(E_long_cheb_k, E_direct_k, atol=1e-8)
end

@testset "energy_long Q2D" begin

    @info "testing the long range part of the thin system"
    n_atoms = 32
    L = (100.0, 100.0, 0.5)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    N_grid = (32, 32, 64)
    uspara = USeriesPara(2)
    M_mid = 0

    k_x, k_y, r_z, H_r, H_c, phase_x, phase_y,  sort_z, z = long_paras_gen(L, N_grid, n_atoms)
    cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, 64)

    E_long_loop_k = energy_long_loop_k(qs, poses, L, M_mid, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, uspara)
    E_long_cheb_k = energy_long_cheb_k(qs, poses, L, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, cheb_mat)

    E_direct_k = long_energy_us_k(qs, poses, 30, L, uspara, M_mid + 1, length(uspara.sw))
    
    @test isapprox(E_long_loop_k, E_long_cheb_k)
    @test isapprox(E_long_loop_k, E_direct_k, atol=1e-8)
    @test isapprox(E_long_cheb_k, E_direct_k, atol=1e-8)
end
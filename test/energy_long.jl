@testset "energy_long long range" begin
    @info "testing the long range part of the energy"
    n_atoms = 8
    L = (100.0, 100.0, 100.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    # using the 3DFFT
    N_grid = (32, 32, 32)
    uspara = USeriesPara(2)
    soepara = SoePara16()
    M_mid = 5

    k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid, USeriesPara(2), M_mid, n_atoms)

    @info "running the FFCT for the long range part of the energy"
    E_FFCT_loop_k = energy_long_loop_k(qs, poses, L, N_grid, M_mid, k_x, k_y, r_z, phase_x, phase_y, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, uspara)
    
    E_FFCT_einsum_k = energy_long_einsum_k(qs, poses, L, M_mid, k_x, k_y, k_mat, r_z, phase_x, phase_y, phase_xs, phase_ys, phase_xys, temp_ijlk, temp_ijl, size_dict, z_coef, exp_coef, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, uspara)

    E_FFCT_0 = energy_long_0(qs, poses, L, M_mid, z, sort_z, uspara, soepara)

    @info "running the direct summation for the long range part of the energy"
    # using the direct summation
    E_direct_k = long_energy_us_k(qs, poses, 30, L, uspara, M_mid + 1, length(uspara.sw))
    E_direct_0 = long_energy_us_0(qs, poses, L, uspara, M_mid + 1, length(uspara.sw))

    @show E_FFCT_loop_k, E_FFCT_einsum_k, E_direct_k
    @test E_FFCT_0 ≈ E_direct_0
    @test E_FFCT_einsum_k ≈ E_FFCT_loop_k
    @test isapprox(E_FFCT_loop_k, E_direct_k, atol=1e-4)
    @test isapprox(E_FFCT_einsum_k, E_direct_k, atol=1e-4)
end
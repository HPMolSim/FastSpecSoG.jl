@testset "interpolate on non-uniform grid" begin
    n_atoms = 8
    L = (100.0, 100.0, 100.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    # using the 3DFFT
    N_grid = (16, 16, 16)
    uspara = USeriesPara(2)
    soepara = SoePara16()
    
    for M_mid in 0:length(uspara.sw) - 1
        @testset "M_mid = $M_mid" begin
            k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid, USeriesPara(2), M_mid, n_atoms)

            H_r_einsum = copy(interpolate_nu_einsum!(copy(H_r), qs, poses, L, k_x, k_y, k_mat, phase_xs, phase_ys, phase_xys, z_coef, exp_coef, temp_ijlk, temp_ijl, r_z, us_mat, uspara, M_mid, size_dict))

            H_r_einsum_direct = copy(interpolate_nu_einsum_non_inplace!(copy(H_r), qs, poses, L, k_x, k_y, k_mat, phase_xs, phase_ys, phase_xys, z_coef, exp_coef, r_z, us_mat, uspara, M_mid))

            H_r_loop = copy(interpolate_nu_loop!(copy(H_r), qs, poses, N_grid, L, k_x, k_y, phase_x, phase_y, r_z, us_mat, uspara, M_mid))

            @test H_r_einsum ≈ H_r_einsum_direct
            @test H_r_einsum_direct ≈ H_r_loop
            @test H_r_einsum ≈ H_r_loop
        end
    end
end
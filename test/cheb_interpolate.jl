@testset "Cheb expansion of U_series" begin
    @info "testing the Cheb expansion"
    k_xi = 0.1
    k_yj = 0.2
    L_z = 10.0
    z = [0.01:0.01:L_z...]
    f_cheb = FastSpecSoG.ChebUSeries_k(k_xi, k_yj, L_z, USeriesPara(2), 3, 16)
    direct_U = map(x -> FastSpecSoG.USeries_direct(x, k_xi, k_yj, USeriesPara(2), 3), z)
    @test isapprox(f_cheb.(z), direct_U)
end

@testset "real2cheb and cheb2real" begin
    n_atoms = 32
    L = (100.0, 100.0, 1.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    N_grid = (32, 32, 32)
    uspara = USeriesPara(2)
    M_mid = 0

    k_x, k_y, k_mat, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, phase_x, phase_y, phase_xs, phase_ys, phase_xys, rhs, sol, sort_z, z, size_dict, temp_ijlk, temp_ijl, z_coef, exp_coef = FFCT_precompute(L, N_grid, uspara, M_mid, n_atoms)

    r_z = chebpoints(N_grid[3], 0.0, L[3])
    H_r = interpolate_Q2D_direct!(H_r, qs, poses, L, k_x, k_y, phase_x, phase_y, r_z, uspara, M_mid)
    H_c = real2Cheb_Q2D!(H_r, H_c, r_z, L[3])
    H_r2 = Cheb2real_Q2D!(H_c, copy(H_r), r_z, L[3])

    @test isapprox(H_r, H_r2; atol = 1e-13)
end
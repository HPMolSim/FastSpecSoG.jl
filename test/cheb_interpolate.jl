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

@testset "real2Cheb and Cheb2real" begin
    n_atoms = 32
    L = (100.0, 100.0, 1.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    N_grid = (32, 32, 32)
    uspara = USeriesPara(2)
    M_mid = 0

    k_x, k_y, r_z, H_r1, H_c, phase_x, phase_y,  sort_z, z = long_paras_gen(L, N_grid, n_atoms)
    cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, 64)

    H_r2 = copy(H_r1)

    H_r = interpolate_nu_loop!(H_r1, qs, poses, L, k_x, k_y, phase_x, phase_y, r_z, uspara, M_mid)
    H_c = real2Cheb!(H_r1, H_c, r_z, L[3])
    H_r2 = Cheb2real!(H_c, H_r1, r_z, L[3])

    @test isapprox(H_r, H_r2; rtol = 1e-14)
end
@testset "interpolate on non-uniform grid" begin
    @info "testing the interpolation on non-uniform grid"
    n_atoms = 8
    L = (100.0, 100.0, 100.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    # using the FFCT
    N_grid = (16, 16, 16)
    uspara = USeriesPara(2)
    Q = 64

    for M_mid in 2:length(uspara.sw) - 1
        @testset "M_mid = $M_mid" begin
            k_x, k_y, r_z, H_r, H_c, phase_x, phase_y = long_paras_gen(L, N_grid)
            cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, Q)

            H_r_loop = copy(interpolate_nu_loop!(copy(H_r), qs, poses, L, k_x, k_y, phase_x, phase_y, r_z, uspara, M_mid))

            H_r_cheb = copy(interpolate_nu_cheb!(copy(H_r), qs, poses, L, k_x, k_y, phase_x, phase_y, r_z, cheb_mat))

            @test isapprox(H_r_loop, H_r_cheb; atol=1e-8)
        end
    end
end
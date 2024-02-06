@testset "energy_long long range" begin
    @info "testing the long range part of the energy"
    n_atoms = 8
    L = (100.0, 100.0, 100.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    # using the 3DFFT
    N_grid = (16, 18, 10)
    uspara = USeriesPara(2)
    soepara = SoePara16()
    M_mid = 8

    k_x, k_y, r_z, us_mat, H_r, H_c, H_s, ivsm, b_l, b_u, rhs, sol, sort_z, z = FFCT_precompute(L, N_grid, USeriesPara(2), M_mid, n_atoms)

    @info "running the FFCT for the long range part of the energy"
    E_FFCT = energy_long(qs, poses, L, M_mid, k_x, k_y, r_z, z, sort_z, us_mat, b_l, b_u, rhs, sol, ivsm, H_r, H_c, H_s, uspara, soepara)

    @info "running the direct summation for the long range part of the energy"
    # using the direct summation
    E_direct = long_energy_us(qs, poses, 30, L, uspara, M_mid + 1, length(uspara.sw))

    @show E_FFCT, E_direct
    @test isapprox(E_FFCT, E_direct, atol=1e-4)
end
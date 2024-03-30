@testset "thin systems long" begin

    @info "testing the long range part of the thin system, 2dfft"

    preset = 2
    
    n_atoms = 32
    L = (100.0, 100.0, 1.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    # 2dfft
    N_real = (128, 128)
    R_z = 16
    w = (16, 16)
    β = 5.0 .* w
    cheb_order = 16
    uspara = USeriesPara(preset)
    Taylor_Q = 16

    gridinfo, pad_grids, cheb_coefs, scalefactors, H_r, H_c, cheb_value, r_z = thin_paras_gen(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q)

    E_thin_k = energy_long_thin_k(qs, poses, L, r_z, H_r, H_c, gridinfo, pad_grids, scalefactors, cheb_coefs, cheb_value)

    E_direct_k = long_energy_us_k(qs, poses, 30, L, uspara, 1, length(uspara.sw))

    @test isapprox(E_thin_k, E_direct_k)
end

@testset "energy_long Q2D" begin

    @info "testing the long range part of the thin system, loop and cheb"
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
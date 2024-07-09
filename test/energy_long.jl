@testset "energy long range" begin
    @info "testing the long range part of the energy"
    n_atoms = 32
    L = (100.0, 100.0, 100.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    N_grid = (32, 32, 32)
    uspara = USeriesPara(2)
    M_mid = 5

    k_x, k_y, r_z, H_r, H_c, phase_x, phase_y = long_paras_gen(L, N_grid)
    cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, 64)
    r_z0, grids0, chebcoefs0 = zero_paras_gen(L[3], 64)
    chebuseies = ChebUSeries_0(L[3], uspara, M_mid, 64)

    E_long_loop_k = energy_long_loop_k(qs, poses, L, M_mid, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, uspara)
    E_long_cheb_k = energy_long_cheb_k(qs, poses, L, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, cheb_mat)

    E_long_0 = energy_long_0(qs, poses, L, chebuseies, r_z0, chebcoefs0, grids0)

    @info "running the direct summation for the long range part of the energy"
    # using the direct summation
    E_direct_k = long_energy_us_k(qs, poses, 10^(-16), L, uspara, M_mid + 1, length(uspara.sw))
    E_direct_0 = long_energy_us_0(qs, poses, L, uspara, M_mid + 1, length(uspara.sw))

    @test isapprox(E_long_0, E_direct_0, atol = 1e-4)
    @test isapprox(E_long_loop_k, E_long_cheb_k)
    @test isapprox(E_long_loop_k, E_direct_k, atol=1e-8)
    @test isapprox(E_long_cheb_k, E_direct_k, atol=1e-8)
end

@testset "thin systems long, fft" begin

    @info "testing the long range part of the thin system, 2dfft"

    preset = 1
    
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

    E_direct_k = long_energy_us_k(qs, poses, 10^(-16), L, uspara, 1, length(uspara.sw))

    @test isapprox(E_thin_k, E_direct_k)
end

@testset "thin systems long, loop & cheb" begin

    @info "testing the long range part of the thin system, loop and cheb"
    n_atoms = 32
    L = (100.0, 100.0, 0.1)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    N_grid = (32, 32, 64)
    uspara = USeriesPara(1)
    M_mid = 3

    k_x, k_y, r_z, H_r, H_c, phase_x, phase_y = long_paras_gen(L, N_grid)
    cheb_mat = ChebUSeries(k_x, k_y, L[3], uspara, M_mid, 64)

    E_long_loop_k = energy_long_loop_k(qs, poses, L, M_mid, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, uspara)
    E_long_cheb_k = energy_long_cheb_k(qs, poses, L, k_x, k_y, r_z, phase_x, phase_y, H_r, H_c, cheb_mat)

    E_direct_k = long_energy_us_k(qs, poses, 10^(-16), L, uspara, M_mid + 1, length(uspara.sw))

    slab_nufft = SlabNUFFT(n_atoms, (64, 64), 8, L, 1e-15, uspara, M_mid)
    E_nufft_slab = nufft_energy_long_k(qs, poses, slab_nufft)
    
    @test isapprox(E_long_loop_k, E_long_cheb_k)
    @test isapprox(E_nufft_slab, E_direct_k)
    @test isapprox(E_long_loop_k, E_direct_k, atol=1e-8)
    @test isapprox(E_long_cheb_k, E_direct_k, atol=1e-8)
end
@testset "thin systems with 2dfft" begin

    preset = 2
    
    n_atoms = 16
    L = (100.0, 100.0, 2.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    N_real = (32, 32)
    R_z = 5
    w = (16, 16)
    β = 5.0 .* w
    cheb_order = 10
    uspara = USeriesPara(preset)
    Taylor_Q = 10

    gridinfo, gridboxs, cheb_coefs, scalefactors, qs_new, poses_new, H_r, H_c, H_s, cheb_value, r_z = thin_precompute(N_real, R_z, w, β, L, cheb_order, uspara, Taylor_Q, n_atoms)

    E_thin_k = energy_thin_k(qs, poses, qs_new, poses_new, L, r_z, H_r, H_c, gridinfo, gridboxs, cheb_coefs, scalefactors, cheb_value)

    E_direct_k = long_energy_us_k(qs, poses, 30, L, uspara, 1, length(uspara.sw))

    @show E_thin_k, E_direct_k
    @test isapprox(E_thin_k, E_direct_k, atol=1e-10)
end
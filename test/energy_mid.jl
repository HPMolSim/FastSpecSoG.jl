@testset "energy_long middle range" begin
    @info "testing the middle range part of the energy"
    n_atoms = 8
    L = (100.0, 100.0, 100.0)

    qs = rand(n_atoms)
    qs .-= sum(qs) ./ n_atoms
    poses = [tuple(L .* rand(3)...) for i in 1:n_atoms]

    # using the 3DFFT
    N_real = (32, 32, 32)
    w = Int.(N_real ./ 4)
    β = 5.0 .* w
    extra_pad_ratio = 2
    cheb_order = 10
    uspara = USeriesPara(2)
    M_mid = 3

    gridinfo, gridbox, cheb_coefs, scalefactor = mid_paras_gen(N_real, w, β, L, extra_pad_ratio, cheb_order, uspara, M_mid)

    @info "running the 3DFFT for the middle range part of the energy"
    E_3DFFT = energy_mid(qs, poses, gridinfo, gridbox, cheb_coefs, scalefactor)

    @info "running the direct summation for the middle range part of the energy"
    # using the direct summation
    E_direct = long_energy_us(qs, poses, 30, L, uspara, 1, M_mid)

    @show E_3DFFT, E_direct
    @test isapprox(E_3DFFT, E_direct, atol=1e-4)
end
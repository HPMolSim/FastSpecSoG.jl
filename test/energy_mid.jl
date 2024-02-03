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
    E_direct = 0.0
    cutoff = 30
    for i in ProgressBar(1:n_atoms)
        for j in 1:n_atoms
            for mx in -cutoff:cutoff
                for my in -cutoff:cutoff
                    for l in 1:M_mid
                        sl, wl = uspara.sw[l]
                        E_direct += qs[i] * qs[j] * π / (L[1] * L[2]) * wl * sl^2 * exp(- (poses[i][3] - poses[j][3])^2 / sl^2) * exp(- sl^2 * ((2π * mx / L[1])^2 + (2π * my / L[2])^2)/4) * cos(2π * mx / L[1] * (poses[i][1] - poses[j][1]) + 2π * my / L[2] * (poses[i][2] - poses[j][2]))
                    end
                end
            end
        end
    end

    @show E_3DFFT, E_direct
    @test isapprox(E_3DFFT, E_direct, atol=1e-4)
end
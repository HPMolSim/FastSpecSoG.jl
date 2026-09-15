@testset "core works without ExTinyMD" begin
    # The requirement this whole phase exists for: FastSpecSoG's plans must be
    # constructible and queryable with ExTinyMD never loaded. A @testset inside
    # the normal suite cannot prove that on its own -- test/runtests.jl itself
    # does `using ExTinyMD` (for test/adapter.jl and for the Ewald2D reference
    # in test/energy.jl), so by the time this runs, ExTinyMD is already loaded
    # in THIS process. The only way to prove the core does not need it is a
    # fresh process that never imports it.
    script = """
    using FastSpecSoG
    @assert !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD"))

    n = 24
    L = (20.0, 20.0, 20.0)
    r_c = 4.0      # r_c = 4.0 < min(Lx, Ly)/2 = 10.0

    # A plain lattice: no RNG needed, every pair separation is an exact
    # multiple of the spacing (including one pair whose minimum image comes
    # from an exact x wrap), which is the configuration most likely to expose a
    # 0/0 in the short-range sum. Two z layers so the out-of-plane part of the
    # separation is exercised.
    poses = [(i * 2.5, j * 3.0, 4.0 + (i % 2) * 2.0) for i in 0:7 for j in 0:2]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    @assert length(poses) == n
    @assert sum(charges) == 0.0

    naive = FSSoG_naive(L, n, r_c, 2.0, preset = 3)
    E_naive = energy_naive(naive, poses, charges)
    E_long_naive = long_energy_naive(naive, poses, charges)
    E_short_naive = short_energy_naive(naive, poses, charges)

    uspara_cheb, F0 = Es_Cheb_precompute(3, 0.5, r_c, 64)
    E_cheb = short_energy_Cheb(uspara_cheb, r_c, F0, L, poses, charges)

    plan = FSSoGInteraction(L, n, r_c, 48, 0.5, (64, 64, 64), (8, 8, 8),
                            5.0 .* (8, 8, 8), 2, 10, 3, (16, 16, 16), 48, 16, 16;
                            preset = 3, ϵ = 1.0)
    E = FastSpecSoG.energy(plan, poses, charges)
    E_pa = energy_per_atom(plan, poses, charges)

    thin = FSSoGThinInteraction((20.0, 20.0, 1.0), n, r_c, 48, 0.5, (64, 64), 16,
                                (8, 8), 5.0 .* (8, 8), 16, 24, 16, 16;
                                preset = 3, ϵ = 1.0)
    poses_thin = [(p[1], p[2], 0.1 + 0.4 * (p[3] - 4.0)) for p in poses]
    E_thin = FastSpecSoG.energy(thin, poses_thin, charges)

    # A candidate-pair list, built by hand: no neighbour-finder type needed.
    nb = [(i, j, 0.0) for i in 1:n-1 for j in i+1:n]
    E_nb = FastSpecSoG.energy(plan, poses, charges; neighbor_list = nb)

    for v in (E_naive, E_long_naive, E_short_naive, E_cheb, E, E_thin, E_nb)
        @assert isfinite(v)
    end
    @assert all(isfinite, E_pa)
    @assert isapprox(E_nb, E; rtol = 1e-14)
    @assert isapprox(sum(E_pa), E; atol = 1e-12)

    # Still not loaded after every query path has run: nothing in src/ pulls it
    # in lazily.
    @assert !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD"))

    print("OK")
    """
    out = read(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script`, String)
    @test out == "OK"

    # Verify the guard above is not vacuous: it must actually fail when
    # ExTinyMD IS loaded first.
    #
    # `@test !success(proc)` alone would NOT be enough, and that is the whole
    # point: it is true for any non-zero exit -- ExTinyMD missing from the test
    # environment, an unsatisfiable resolve, a failed precompile, a typo in the
    # heredoc. Every one of those makes the check pass while proving nothing
    # about whether the `@assert` on `Base.loaded_modules` is load-bearing.
    # (That was a review finding on QuasiEwald in this same phase.) So the
    # subprocess's streams are captured and the failure is pinned to the
    # intended cause: a marker printed after the `using` lines must appear on
    # stdout (so the loads all succeeded), the marker after the assert must NOT
    # appear (so the assert is what stopped it), and stderr must carry an
    # AssertionError naming the guard's own expression.
    script_contaminated = """
    using ExTinyMD, FastSpecSoG
    print("LOADS_OK;")
    @assert !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD"))
    print("ASSERT_DID_NOT_TRIP;")
    """
    outfile, errfile = tempname(), tempname()
    proc = run(pipeline(
        `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script_contaminated`;
        stdout = outfile, stderr = errfile), wait = false)
    wait(proc)
    out_text = read(outfile, String)
    err_text = read(errfile, String)
    rm(outfile, force = true)
    rm(errfile, force = true)

    @test !success(proc)
    # `using ExTinyMD, FastSpecSoG` really did succeed, so the non-zero exit is
    # not a missing package, a resolve failure or a precompilation error.
    @test occursin("LOADS_OK;", out_text)
    # and execution stopped at the assert, not after it.
    @test !occursin("ASSERT_DID_NOT_TRIP;", out_text)
    # and it stopped for exactly the intended reason.
    @test occursin("AssertionError", err_text)
    @test occursin("loaded_modules", err_text)
end

# Scan a directory tree of .jl files for `ExTinyMD` appearing in CODE. Comments
# and docstrings may mention it -- several in src/ explain what was replaced and
# why -- so comment-only lines, trailing comments and the bodies of `"""..."""`
# docstrings are stripped first.
function _extinymd_code_refs(dir)
    offenders = String[]
    for (root, _, files) in walkdir(dir), f in files
        endswith(f, ".jl") || continue
        path = joinpath(root, f)
        in_docstring = false
        for (lineno, line) in enumerate(eachline(path))
            n_quotes = count("\"\"\"", line)
            was_in = in_docstring
            isodd(n_quotes) && (in_docstring = !in_docstring)
            (was_in || in_docstring) && continue
            stripped = strip(line)
            startswith(stripped, "#") && continue
            code = split(line, '#')[1]
            occursin("ExTinyMD", code) && push!(offenders, "$path:$lineno: $stripped")
        end
    end
    return offenders
end

@testset "src/ has no ExTinyMD code references" begin
    # The standing constraint, checked mechanically rather than by review.
    root = dirname(@__DIR__)
    offenders = _extinymd_code_refs(joinpath(root, "src"))
    isempty(offenders) || @info "ExTinyMD references found in src/ code" offenders
    @test isempty(offenders)

    # The scanner is not vacuous: run on ext/, which SHOULD be full of them, it
    # finds them. Without this, a scanner that silently matched nothing would
    # pass the assertion above forever.
    @test !isempty(_extinymd_code_refs(joinpath(root, "ext")))

    # ExTinyMD is a weak dependency only.
    project = read(joinpath(root, "Project.toml"), String)
    deps = split(split(project, "[deps]")[2], "[weakdeps]")[1]
    @test !occursin("ExTinyMD", deps)
    weak = split(split(project, "[weakdeps]")[2], "[extensions]")[1]
    @test occursin("ExTinyMD", weak)

    # Every r_c literal in test/ carries an inline "< min(Lx, Ly)/2 = ..."
    # annotation, which is what documents the minimum-image requirement at each
    # site. Counted rather than eyeballed, so removing one is visible. 15 sites
    # at the time of writing; the bound is 14 so that adding a test file does
    # not require touching this number, while deleting two annotations does.
    n_annotated = 0
    for f in readdir(joinpath(root, "test"); join = true)
        endswith(f, ".jl") || continue
        for line in eachline(f)
            # skip this scanner's own source lines
            occursin("occursin(", line) && continue
            occursin("min(Lx, Ly)/2", line) || continue
            n_annotated += 1
        end
    end
    @test n_annotated >= 14
end

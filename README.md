# FastSpecSoG.jl

[![Build Status](https://github.com/HPMolSim/FastSpecSoG.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HPMolSim/FastSpecSoG.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HPMolSim/FastSpecSoG.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HPMolSim/FastSpecSoG.jl)


`FastSpecSoG.jl` is an implementation of the newly developed fast spectral sum-of-Gaussian method for the calculation of the electrostatic potential and field in double periodic molecular systems in Julia Programming Language. 
The method is based on the sum-of-Gaussian (SoG) approximation, the Non-Uniform Fast Fourier Transform (NUFFT) algorithm and the FFT-FFCT mixed algorithm and reached a complexity of $O(N \log N)$ with spectral accuracy.

For more details about our method, please refer to our arxiv article [A fast spectral sum-of-Gaussians method for electrostatic summation in quasi-2D systems]([https://arxiv.org/abs/2412.04595](https://doi.org/10.1007/s00211-025-01518-y)).
For benchmarks, please see this repo https://github.com/xuanzhaogao/FSSoG_benchmark.

## Getting Started

```julia
pkg> add FastSpecSoG
```

`ExTinyMD` is **not** required. It is a weak dependency: install it only if you
want the MD adapter, and load it before you use it.

```julia
pkg> add ExTinyMD        # optional
```

## Standalone usage (no ExTinyMD)

Everything in `src/` is framework-free. A plan is built from the box and the
particle count, and then queried with plain arrays:

```julia
using FastSpecSoG

n_atoms = 100
L = (50.0, 50.0, 50.0)
r_c = 10.0      # must be strictly less than min(Lx, Ly)/2 = 25.0

# array-of-structs positions and a charge vector -- your own arrays, untouched
poses   = [(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n_atoms]
charges = [isodd(i) ? 2.0 : -2.0 for i in 1:n_atoms]

plan = FSSoGInteraction(L, n_atoms, r_c, 48, 0.5, (128, 128, 128), (16, 16, 16),
                        5.0 .* (16, 16, 16), 2, 10, 3, (32, 32, 32), 48, 32, 32;
                        preset = 3, ϵ = 1.0)

E  = FastSpecSoG.energy(plan, poses, charges)
Ei = energy_per_atom(plan, poses, charges)
```

Four things are worth knowing about this API.

**`poses` is array-of-structs and its element type is not constrained.** The
kernels only ever index `p[1]`, `p[2]`, `p[3]`, so `NTuple{3,T}`,
`SVector{3,T}` and an MD framework's own point type all work directly, with no
conversion layer.

**Your arrays are never modified.** A query copies positions and charges into
scratch the plan owns, so nothing is allocated per call and nothing you passed
in is touched.

**`r_c` must be strictly less than `min(Lx, Ly)/2`.** The short-range sum uses
the minimum image under the x/y periodicity (z is the free slab axis), and the
nearest image is only unique below that bound.

**`energy` is defined but NOT exported.** Several sibling electrostatics
packages define a function of the same name, so exporting it from all of them
would make the bare name ambiguous under `using FastSpecSoG, SomeOther`. Call
it as `FastSpecSoG.energy(...)`. The descriptive names stay exported and are
unchanged: `energy_naive`, `energy_short`, `energy_mid`, `energy_long`,
`energy_per_atom`, `short_energy_naive`, `short_energy_Cheb`,
`long_energy_naive`, and the rest.

### Passing a neighbour list

The short-range sum takes all `i < j` pairs by default, which is `O(N^2)`. If
you already maintain a neighbour list, pass it:

```julia
E = FastSpecSoG.energy(plan, poses, charges; neighbor_list = my_list)
```

`my_list` is any iterable whose elements support `pair[1]` and `pair[2]`, or an
object that owns such a list under a `neighbor_list` property (which is how MD
cell lists usually expose theirs). **Only the indices are used.** Any distance
the list carries is ignored and the true three-dimensional separation is
recomputed from the positions, because a quasi-2D cell list reports the
in-plane distance, not the separation the energy needs. An in-plane list is a
superset of the true pair set, so re-filtering on the recomputed distance gives
the correct answer rather than merely rejecting bad data.

Pairs at exactly zero separation -- two particles at one site, or two separated
by an exact multiple of `Lx` or `Ly` with matching `y` and `z`, as a lattice
initialisation produces -- are skipped, since the short-range kernel divides by
`r`.

## MD usage via ExTinyMD

Loading `ExTinyMD` alongside this package activates
`ext/FastSpecSoGExTinyMDExt.jl`, which adds one method:

```julia
using FastSpecSoG, ExTinyMD

E  = ExTinyMD.energy(plan, neighborfinder, sys, info)
Ei = FastSpecSoG.energy_per_atom(plan, neighborfinder, sys, info)
```

It gathers positions and charges in storage-slot order, honouring ExTinyMD's
id/slot indirection (`sys.atoms` is indexed by particle id,
`info.particle_info` by storage slot), and calls the core.

**FastSpecSoG is energy-only and cannot drive an MD run.** It has never
computed forces: there is no `force`, no `force!` and no
`update_acceleration!`. `simulate!` calls `update_acceleration!` on every
interaction at every step, so a FastSpecSoG plan could not be integrated
regardless of this adapter.

Relatedly, the plan types cannot be placed in `sys.interactions` or in an
`EnergyLogger`: both require `ExTinyMD.AbstractInteraction`, and a struct's
supertype is fixed where the struct is defined -- these are defined in `src/`,
which has no ExTinyMD, so no extension can retrofit one. `ExTinyMD.energy` on a
FastSpecSoG plan is therefore callable **directly and only directly**. For
electrostatics that can drive `simulate!`, use ExTinyMD's own `Ewald2D` or
`PME3D`, or `QuasiEwald.jl`/`SoEwald2D.jl`.

### Changed from 0.1.x

`ExTinyMD` moved from `[deps]` to `[weakdeps]`, and the query API changed, so
0.2.0 is a breaking release.

| 0.1.x | 0.2.0 |
|---|---|
| `ExTinyMD.energy(plan, neighbor, info, atoms)` | `ExTinyMD.energy(plan, finder, sys, info)`, or `FastSpecSoG.energy(plan, poses, charges)` |
| `energy_naive(plan, neighbor, info, atoms)` | `energy_naive(plan, poses, charges; neighbor_list)` |
| `energy_per_atom(plan, neighbor, info, atoms)` | `energy_per_atom(plan, poses, charges; neighbor_list)` |
| `short_energy_naive(plan, neighbor, position, q)` | `short_energy_naive(plan, position, q; neighbor_list)` |
| `short_energy_Cheb(cheb, r_c, F0, boundary, neighbor, position, q)` | `short_energy_Cheb(cheb, r_c, F0, L, position, q; neighbor_list)` |
| `energy_short(plan, neighbor)` | `energy_short(plan; neighbor_list)` |
| `plan.boundary` | removed -- it was only ever `Q2dBoundary(plan.L...)` |

The old `ExTinyMD.energy(plan, neighbor, info, atoms)` methods were already
unreachable from ExTinyMD: its MD loop calls
`energy(interaction, neighborfinder, sys, info)`, a different argument order,
so nothing but this package's own tests ever called them. The replacement
carries the signature ExTinyMD actually uses.

## Numerical Methods

In this package, two sets of method are implemented. Consider a simulation box with edge length $L_x$, $L_y$ and $L_z$ and $N$ atoms, with is double periodic in $xy$ and non-periodic in $z$ direction.
We assume $L_x \approx L_y$, and if $L_z \approx L_x$, we call the system a cubic system, and if $L_z \ll L_x$, we call the system a slab system.

For cubic systems, please refer to the following example:
```julia
using FastSpecSoG

n_atoms = 100
L = (50.0, 50.0, 50.0)

poses   = [(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n_atoms]
charges = [isodd(i) ? 2.0 : -2.0 for i in 1:n_atoms]

r_c = 10.0      # < min(Lx, Ly)/2 = 25.0
N_real = (128, 128, 128)
w = (16, 16, 16)
β = 5.0 .* w
extra_pad_ratio = 2
cheb_order = 10
preset = 3
M_mid = 3

N_grid = (32, 32, 32)
Q = 48
R_z0 = 32
Q_0 = 32

fssog_interaction = FSSoGInteraction(L, n_atoms, r_c, Q, 0.5, N_real, w, β,
                                     extra_pad_ratio, cheb_order, M_mid,
                                     N_grid, Q, R_z0, Q_0;
                                     preset = preset, ϵ = 1.0)

energy_sog = FastSpecSoG.energy(fssog_interaction, poses, charges)
```

For slab systems, please refer to the following example
```julia
using FastSpecSoG

n_atoms = 100
L = (100.0, 100.0, 1.0)

poses   = [(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n_atoms]
charges = [isodd(i) ? 2.0 : -2.0 for i in 1:n_atoms]

r_c = 10.0      # < min(Lx, Ly)/2 = 50.0
N_real = (128, 128)
R_z = 32
w = (16, 16)
β = 5.0 .* w
cheb_order = 16
preset = 3
Q = 48
Q_0 = 32
R_z0 = 32
Taylor_Q = 24

fssog_interaction = FSSoGThinInteraction(L, n_atoms, r_c, Q, 0.5, N_real, R_z,
                                         w, β, cheb_order, Taylor_Q, R_z0, Q_0;
                                         preset = preset, ϵ = 1.0)

energy_sog = FastSpecSoG.energy(fssog_interaction, poses, charges)
```

Both examples place particles at random, which can put two closer than the
`r_min = 0.5` floor of the short-range Chebyshev interpolant; use a
configuration with a minimum separation (or raise `r_min`) for real work.

## Citation

If you use this package in your research, please cite our article:
```
@article{fssog,
      title={{A fast spectral sum-of-Gaussians method for electrostatic summation in quasi-2D systems}}, 
      author={Xuanzhao Gao and Shidong Jiang and Jiuyang Liang and Zhenli Xu and Qi Zhou},
      year={2024},
      journal={arXiv preprint arXiv:2412.04595},
}
```

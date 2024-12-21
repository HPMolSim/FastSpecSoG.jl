# FastSpecSoG.jl

[![Build Status](https://github.com/HPMolSim/FastSpecSoG.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HPMolSim/FastSpecSoG.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HPMolSim/FastSpecSoG.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HPMolSim/FastSpecSoG.jl)


`FastSpecSoG.jl` is an implementation of the newly developed fast spectral sum-of-Gaussian method for the calculation of the electrostatic potential and field in double periodic molecular systems in Julia Programming Language. 
The method is based on the sum-of-Gaussian (SoG) approximation, the Non-Uniform Fast Fourier Transform (NUFFT) algorithm and the FFT-FFCT mixed algorithm and reached a complexity of $O(N \log N)$ with spectral accuracy.

For more details about our method, please refer to our arxiv article [A fast spectral sum-of-Gaussians method for electrostatic summation in quasi-2D systems](https://arxiv.org/abs/2412.04595).
For benchmarks, please see this repo https://github.com/ArrogantGao/FSSoG_benchmark.

## Getting Started

First you need to add its deps `ExTinyMD` in `Julia` by typing `]` in Julia REPL and then
```julia
pkg> add ExTinyMD
```
Then you can add this package by typing
```julia
pkg> add FastSpecSoG
```
to install the package.

## Numerical Methods

In this package, two sets of method are implemented. Consider a simulation box with edge length $L_x$, $L_y$ and $L_z$ and $N$ atoms, with is double periodic in $xy$ and non-periodic in $z$ direction.
We assume $L_x \approx L_y$, and if $L_z \approx L_x$, we call the system a cubic system, and if $L_z \ll L_x$, we call the system a slab system.

For cubic systems, please refer to the following example:
```julia
n_atoms = 100
L = 50.0
boundary = ExTinyMD.Q2dBoundary(L, L, L)

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms÷2
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 2.0))
end

for i in n_atoms÷2 + 1 : n_atoms
    push!(atoms, Atom(type = 2, mass = 1.0, charge = - 2.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

r_c = 10.0
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

fssog_interaction = FSSoGInteraction((L, L, L), n_atoms, r_c, Q, 0.5, N_real, w, β, extra_pad_ratio, cheb_order, M_mid, N_grid, Q, R_z0, Q_0; preset = preset, ϵ = 1.0)

fssog_neighbor = CellList3D(info, fssog_interaction.r_c, fssog_interaction.boundary, 1)
energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)
```

For slab systems, please refer to the following example
```julia
n_atoms = 100
Lx = 100.0
Ly = 100.0
Lz = 1.0
boundary = ExTinyMD.Q2dBoundary(Lx, Ly, Lz)

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms÷2
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 2.0))
end

for i in n_atoms÷2 + 1 : n_atoms
    push!(atoms, Atom(type = 2, mass = 1.0, charge = - 2.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, Lx, 0.0, Ly, 0.0, Lz), boundary; min_r = 1.0, temp = 1.0)

r_c = 1.0
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

fssog_interaction = FSSoGThinInteraction((Lx, Ly, Lz), n_atoms, r_c, Q, 0.5, N_real, R_z, w, β, cheb_order, Taylor_Q, R_z0, Q_0; preset = preset, ϵ = 1.0)

fssog_neighbor = CellList3D(info, fssog_interaction.r_c, fssog_interaction.boundary, 1)
energy_sog = ExTinyMD.energy(fssog_interaction, fssog_neighbor, info, atoms)
```

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

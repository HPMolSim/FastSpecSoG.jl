# N_x, N_y are the uniform grid numbers in the x and y directions
# N_z is the Chebyshev grid number in the z direction
# n_atoms is the number of particles
# M_num is the number of Gaussians M_mid+1:M

# This part will be done in the following steps:
# 1. Generate the inverse matrix, k matrix
# 2. interpolate the particles onto the grid (Fourier grids in xy and chebgrid in z)
# 3. Convert the grid in z into the Chebyshev series
# 4. solve the equation for the boundary condition
# 5. gather the energy
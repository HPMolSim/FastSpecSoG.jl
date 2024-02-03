# precomputings

# precompute the parameters to be used
function FFCT_precompute()

    return k_x, k_y, r_z, k_mat, kx_mat, ky_mat
end

# precompute the inverse matrix
function inverse_mat(N_grid::NTuple{3, Int}, L_z::T, k_x::Vector{T}, k_y::Vector{T}) where{T}
    N_x, N_y, N_z = N_grid
    inverse_mat = zeros(T, N_z + 2, N_z + 2, 2 * N_x + 1, 2 * N_y + 1)

    scale = L_z / T(2)

    A = zeros(T, N_z, N_z + 2)
    A[2,2] = 1 / 8 * scale^2
    A[2,4] = - 1 / 8 * scale^2
    for s in 3:N_z
        k = s - 1
        A[s, s - 2] = scale^2 / (2 * k * 2(k - 1))
        A[s, s] = - (1 / 2(k - 1) + 1 / 2(k + 1)) * scale^2 / 2k
        A[s, s + 2] = scale^2 / (2k * 2(k+1))
    end

    B = zeros(T, N_z, 2)
    B[1, 1] = T(2)
    B[2, 2] = scale

    C = zeros(T, 2, N_z + 2)
    for s = 3:N_z
        k = s - 1
        C[1, s-2] += scale^2 * (1/(2*k)) * (1/(2*(k-1))) * (-1)^k
        C[1, s] -= scale^2 * ((1/(2*k)) * (1/(2*(k-1))) + (1/(2*k)) * (1/(2*(k+1)))) * (-1)^k
        C[1, s+2] += scale^2 * (1/(2*k)) * (1/(2*(k+1))) * (-1)^k
        
        C[2, s-2] += scale^2 * (1/(2*k)) * (1/(2*(k-1)))
        C[2, s] -= scale^2 * ((1/(2*k)) * (1/(2*(k-1))) + (1/(2*k)) * (1/(2*(k+1))))
        C[2, s+2] += scale^2 * (1/(2*k)) * (1/(2*(k+1)))
    end
    C[1, 2] += 1/8 * scale^2 * (-1)
    C[1, 4] -= 1/8 * scale^2 * (-1)
    C[2, 2] += 1/8 * scale^2
    C[2, 4] -= 1/8 * scale^2

    D = zeros(T, 2, 2)
    D[1, 1] = T(1)
    D[1, 2] = - scale
    D[2, 1] = T(1)
    D[2, 2] = scale

    divisor = zeros(T, N_z + 2, N_z + 2)
    for i in 1:2N_x + 1, j in 1:2N_y + 1
        k = sqrt(k_x[i]^2 + k_y[j]^2)
        if k != zero(T)
            nu = - k^2
            divisor[1:N_z, 1:N_z] = A[1:N_z, 1:N_z] .* nu .+ T.(I(N_z))
            divisor[1:N_z, N_z+1:N_z+2] = nu .* B[1:N_z, 1:2]
            divisor[N_z+1:N_z+2, 1:N_z] = C[1:2, 1:N_z]
            divisor[N_z+1:N_z+2, N_z+1:N_z+2] = D

            inverse_mat[:,:,i,j] = inv(divisor)
        end
    end

    return inverse_mat
end
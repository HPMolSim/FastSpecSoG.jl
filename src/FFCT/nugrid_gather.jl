@inbounds function gather_nu_single(q::T, pos::NTuple{3, T}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, H_s::Array{Complex{T}, 3}) where{T}
    x, y, z = pos
    L_x, L_y, L_z = L
    ϕ = zero(Complex{T})
    
    revise_phase_pos!(phase_x, phase_y, k_x, k_y, x - T(L[1] / 2), y - T(L[2] / 2))

    for k in 1:size(H_s, 3)
        cheb_val = k == 1 ? T(0.5) : chebpoly(k - 1, z - L_z / T(2), L_z / T(2))
        ϕ += dot(phase_x', (@view H_s[:, :, k]), phase_y)* cheb_val
    end

    return q * ϕ / (2 * L_x * L_y)
end

function gather_nu(qs::Vector{T}, poses::Vector{NTuple{3, T}}, L::NTuple{3, T}, k_x::Vector{T}, k_y::Vector{T}, phase_x::Vector{Complex{T}}, phase_y::Vector{Complex{T}}, H_s::Array{Complex{T}, 3}) where{T}
    E = zero(Complex{T})
    for i in 1:length(qs)
        E += gather_nu_single(qs[i], poses[i], L, k_x, k_y, phase_x, phase_y, H_s)
    end
    return real(E)
end
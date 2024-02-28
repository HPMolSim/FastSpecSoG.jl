# taylor expansion

@inline @inbounds function ChebParticleMesh.horner(x::T, taylorseries::TaylorSeriesPara{N, T}) where {N, T}
    res = taylorseries.coefs[N]
    for i in N - 1:-1:1
        res = res * x + taylorseries.coefs[i]
    end
    return res
end

function TaylorSeriesPara(coefs::Vector{T}) where{T}
    return TaylorSeriesPara{length(coefs), T}(tuple(coefs...))
end

function TaylorUSeries(k_x::Vector{T}, k_y::Vector{T}, L_z::T, uspara::USeriesPara{T}, M_mid::Int, Q::Int) where{T}

    taylor_mat = Array{TaylorSeriesPara{Q, T}, 2}(undef, (length(k_x), length(k_y)))

    for i in 1:length(k_x)
        for j in 1:length(k_y)
            coefs = zeros(T, Q)
            k2 = k_x[i]^2 + k_y[j]^2
            for l in M_mid + 1:length(uspara.sw)
                sl, wl = uspara.sw[l]
                for n in 0:Q - 2
                    coefs[n + 1] += wl / (4π * factorial(n)) * (one(T) + (4n + k2) * sl^2 / 2) * exp(-sl^2 * k2 / 4) * (L_z^2 / sl^2)^n
                end
                coefs[Q] += wl * sl^2 / (2π * factorial(Q - 2)) * exp(-sl^2 * k2 / 4) * (L_z^2 / sl^2)^(Q - 1)
            end
            taylor_mat[i, j] = TaylorSeriesPara(coefs)
        end
    end

    return taylor_mat
end
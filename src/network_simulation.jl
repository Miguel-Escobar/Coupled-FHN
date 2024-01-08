using DifferentialEquations
using StaticArrays
using Statistics

function fhn_eom(x, params)
    a = params[1]
    eps = params[2]
    dx = (x[1] - (x[1]^3)/3 - x[2])/eps
    dy = x[1] + a
    return SVector(dx, dy)
end

function bmatrix(phi, eps)
    return [cos(phi)/eps sin(phi)/eps; -sin(phi) cos(phi)]
end

function coupled_fhn_eom!(dx, x, a, eps, coupling_strength, coupling_matrix, coupling_jac)
    N = length(coupling_matrix[1, :])
    eachneuron = reshape(x, (2, N))
    coupling_terms = coupling_jac * eachneuron
    for i in range(1, N)
        dx_i = fhn_eom(eachneuron[:, i], [a, eps]) .+ coupling_strength .* sum([coupling_matrix[i, j]  .* coupling_terms[:, j] for j in 1:N])
        dx[2*i-1:2*i] = dx_i
    end
    nothing
end

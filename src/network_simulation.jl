using DifferentialEquations
using StaticArrays
using Statistics

include("network_matrices_creation.jl")
include("network_sol_analysis.jl")

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

function coupled_fhn_eom!(dx, x, a, eps, coupling_strength, coupling_matrix, coupling_jac) # not controlled
    N = length(coupling_matrix[1, :])
    eachneuron = reshape(x, (2, N))
    coupling_terms = coupling_jac * eachneuron
    for i in range(1, N)
        dx_i = fhn_eom(eachneuron[:, i], [a, eps]) .+ coupling_strength .* sum([coupling_matrix[i, j]  .* coupling_terms[:, j] for j in 1:N])
        dx[2*i-1:2*i] = dx_i
    end
end


N = 9
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)
σ = 0.0506
G = wattsstrogatzmatrix(N, 3, 0.232)
x_0 = zeros(2*N)
x_0[1:2] .+= 0.1
prob = ODEProblem((dx, x, params, t) -> coupled_fhn_eom!(dx, x, params[1], params[2], params[3], G, b), x_0, (0.0, 200.0), [a, eps, σ])
alg = Tsit5()
sol = solve(prob, alg)
using Plots
plot(sol, xlabel="Time", ylabel="System Variables", dpi=600)
t_val, kuramoto_val = kuramoto_time_series(sol, N)
plot(t_val, kuramoto_val)
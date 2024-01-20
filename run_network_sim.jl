include("src/network_matrices_creation.jl")
include("src/network_simulation.jl")
include("src/network_sol_analysis.jl")
using GLMakie
using LaTeXStrings
using Random
using BenchmarkTools

N = 90
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)
σ = 0.0506
G = wattsstrogatzmatrix(N, 3, 1);#0.232)

x_0 = zeros(2*N)
x_0[2 .* (1:N) .- 1] = rand(N) .* 2 .* a .- a
x_0[2 .* (1:N)] = rand(N) .* 2 .* (-a + a^3 / 3) .- (-a + a^3 / 3)

prob = ODEProblem((dx, x, params, t) -> coupled_fhn_eom!(dx, x, params[1], params[2], params[3], G, b), x_0, (0.0, 500.0), [a, eps, σ])
sol = solve(prob);

@profview solve(prob);
@btime solve(prob);
nothing

# f = Figure(size = (800, 600))
# ax = Axis(f[1, 1])
# ax.title = "Kuramoto Order Parameter"
# ax.xlabel = "Time"
# ax.ylabel = "Kuramoto Order Parameter"
# t_val, kuramoto_val = kuramoto_time_series(sol, N)
# lines!(ax, t_val, kuramoto_val)
# f
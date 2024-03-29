include("src/network_matrices_creation.jl")
include("src/network_simulation.jl")
include("src/network_sol_analysis.jl")
using StatProfilerHTML
using BenchmarkTools
N = 90
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)
σ = 0.0506
G = wattsstrogatzmatrix(N, 3, 1) .* 1.0;#0.232)
# G = convert(Matrix{Float64}, G)
x_0 = zeros(2*N)
x_0[2 .* (1:N) .- 1] = rand(N) .* 2 .* a .- a
x_0[2 .* (1:N)] = rand(N) .* 2 .* (-a + a^3 / 3) .- (-a + a^3 / 3)

# cache = zeros(Float64, N, 2)

dx = zeros(2*N)
cache = zeros(Float64, N, 2)
coupled_fhn_eom!(dx, x_0, a, eps, σ, G, b, N)
x_0[2 .* (1:N) .- 1] = rand(N) .* 2 .* a .- a
x_0[2 .* (1:N)] = rand(N) .* 2 .* (-a + a^3 / 3) .- (-a + a^3 / 3)
@profview for i in 1:100000 coupled_fhn_eom!(dx, x_0, a, eps, σ, G, b, N) end
# @profilehtml for i in 1:10000 coupled_fhn_eom!(dx, x_0, a, eps, σ, G, b) end
@benchmark coupled_fhn_eom!(dx, x_0, a, eps, σ, G, b, N)
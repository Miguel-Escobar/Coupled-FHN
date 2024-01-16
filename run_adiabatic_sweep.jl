include("src/network_matrices_creation.jl")
include("src/cluster_synch.jl")
include("src/network_simulation.jl")
include("src/network_sol_analysis.jl")
include("src/msf.jl")

using Random
using Trapz
using GLMakie
using ProgressBars
using BenchmarkTools
using Base.Threads

N = 10
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)
G = test_matrix_for_cluster_synch()

function kuramoto_d_step(x_0, σ; t_final=1000.0)
    prob = ODEProblem((dx, x, params, t) -> coupled_fhn_eom!(dx, x, params[1], params[2], params[3], G, b), x_0, (0.0, t_final), [a, eps, σ])
    sol = solve(prob; dtmax=0.5)
    t_val, kuramoto_val = kuramoto_time_series(sol, N)
    #x_0 = sol.u[end]
    return trapz(t_val, kuramoto_val)/t_final, sol.u[end]
end

function kuramoto_sweep(start, stop, init_x0, N_d)
    forward_d_sweep = range(start, stop, length=N_d)
    forward_kuramoto_d_vals = zeros(N_d)
    x_0 = init_x0
    for i in 1:N_d
        forward_kuramoto_d_vals[i], x_f = kuramoto_d_step(x_0, forward_d_sweep[i])
        x_0 = x_f
    end
    return forward_d_sweep, forward_kuramoto_d_vals
end


# x_0 = zeros(2*N) .+ randn(2*N) .* 0.01
# kuramoto_d_step!(x_0, 0.015)
# x_0 = zeros(2*N) .+ randn(2*N) .* 0.01
# @btime kuramoto_d_step!(x_0, 0.015)

zero_msf = msf_zero()
println("zero_msf = ", zero_msf)
println("msf = ", master_stability_function(zero_msf, 0))
eigenvalues, eigenvectors, clusters, s_matrices = s_matrix_method(G)
eigenvalues = unique(round.(eigenvalues, digits=8))[2:end]
critical_couplings = zero_msf./eigenvalues #unique(zero_msf./eigenvalues)[2:end]#
println("Critical Couplings: ", critical_couplings)
println("Eigenvalues: ", eigenvalues)

N_d = 100
N_realizations = 100
forward_avg = zeros(N_d)
backward_avg = zeros(N_d)

for i in ProgressBar(1:N_realizations)
    global forward_d_sweep, forward_kuramoto_d = kuramoto_sweep(0.025, 0.005, zeros(2*N) .+ randn(2*N) .* 0.01, N_d)
    global backward_d_sweep, backward_kuramoto_d = kuramoto_sweep(0.005, 0.025, zeros(2*N) .+ randn(2*N) .* 0.01, N_d)
    forward_avg .+= forward_kuramoto_d
    backward_avg .+= backward_kuramoto_d
end
forward_avg ./= N_realizations
backward_avg ./= N_realizations

# backward_d_sweep = reverse(forward_d_sweep)
# backward_kuramoto_d_vals = zeros(N_d)
# init_x_0 = zeros(2*N) .+ randn(2*N) .* 0.01
# for i in ProgressBar(1:N_d)
#     backward_kuramoto_d_vals[i], x_0 = kuramoto_d_step!(x_0, backward_d_sweep[i])
# end

f = Figure(size = (800, 600))
ax = Axis(f[1, 1])
ax.xlabel = "Coupling"
ax.ylabel = "Kuramoto Order Parameter"
scatter!(ax, forward_d_sweep, forward_avg, label="Right to left")
scatter!(ax, backward_d_sweep, backward_avg, label="Left to right")
axislegend()
vlines!(ax, critical_couplings[2:3]; label="Critical Couplings", linewidth=1, color = :red)
f



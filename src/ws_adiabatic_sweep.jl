include("network_matrices_creation.jl")
include("cluster_synch.jl")
include("network_simulation.jl")
include("network_sol_analysis.jl")
include("msf.jl")

using Random
using Trapz
using GLMakie
using ProgressMeter
using BenchmarkTools
using Base.Threads

function kuramoto_d_step(x_0, σ, N, eps, a, b, G; t_final=1000.0)
    prob = ODEProblem((dx, x, params, t) -> coupled_fhn_eom!(dx, x, params[1], params[2], params[3], G, b, N), x_0, (0.0, t_final), [a, eps, σ])
    sol = solve(prob; dtmax=0.5)
    t_val, kuramoto_val = kuramoto_time_series(sol, N)
    #x_0 = sol.u[end]
    return trapz(t_val, kuramoto_val)/t_final, sol.u[end]
end

function kuramoto_sweep(start, stop, init_x0, N_d, N, eps, a, b, G)
    forward_d_sweep = range(start, stop, length=N_d)
    forward_kuramoto_d_vals = zeros(N_d)
    x_0 = init_x0
    @showprogress for i in 1:N_d
        forward_kuramoto_d_vals[i], x_f = kuramoto_d_step(x_0, forward_d_sweep[i], N, eps, a, b, G)
        x_0 = x_f
    end
    return forward_d_sweep, forward_kuramoto_d_vals
end

N = 90
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)

# zero_msf = msf_zero()
# println("zero_msf = ", zero_msf)
# println("msf = ", master_stability_function(zero_msf, 0))
# eigenvalues, eigenvectors, clusters, s_matrices = s_matrix_method(G)
# eigenvalues = unique(round.(eigenvalues, digits=8))[2:end]
# critical_couplings = zero_msf./eigenvalues #unique(zero_msf./eigenvalues)[2:end]#
# println("Critical Couplings: ", critical_couplings)
# println("Eigenvalues: ", eigenvalues)


N_d = 150
N_realizations = 30
forward_array = zeros(N_realizations, N_d)
backward_array = zeros(N_realizations, N_d)

for i in 1:N_realizations
    local G = wattsstrogatzmatrix(N, 3, 1) #test_matrix_for_cluster_synch()
    println("Realization ", i)
    global forward_d_sweep, forward_kuramoto_d = kuramoto_sweep(0.25, 0.00, zeros(2*N) .+ randn(2*N) .* 0.01, N_d, N, eps, a, b, G)
    global backward_d_sweep, backward_kuramoto_d = kuramoto_sweep(0.00, 0.25, zeros(2*N) .+ randn(2*N) .* 0.01, N_d, N, eps, a, b, G)
    forward_array[i, :] .= forward_kuramoto_d
    backward_array[i, :] .= backward_kuramoto_d
end

forward_avg = mean(forward_array, dims=1)[1, :]
backward_avg = mean(backward_array, dims=1)[1, :]

using Serialization
serialize("ws_sweep_data/forward_array", forward_array)
serialize("ws_sweep_data/backward_array", backward_array)
serialize("ws_sweep_data/forward_d_sweep", forward_d_sweep)
serialize("ws_sweep_data/backward_d_sweep", backward_d_sweep) # This will do until I write the poster.


f = Figure(size = (800, 600))
ax = Axis(f[1, 1])
ax.xlabel = "Coupling"
ax.ylabel = "Kuramoto Order Parameter"
scatter!(ax, forward_d_sweep[5:end-5], forward_avg[5:end-5], label="Right to left") # 5000 timesteps for thermalization
scatter!(ax, backward_d_sweep[5:end-5], backward_avg[5:end-5], label="Left to right")
axislegend(position=:rb)
# vlines!(ax, critical_couplings[2:3]; label="Critical Couplings", linewidth=1, color = :red)
# save("coupling_sweep_watts_strogatz_90_neurons_6_neighbours_reconnectionprob_1.png", f)
f

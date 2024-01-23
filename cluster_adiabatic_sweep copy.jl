include("src/network_matrices_creation.jl")
include("src/cluster_synch.jl")
include("src/network_simulation.jl")
include("src/network_sol_analysis.jl")
include("src/msf.jl")

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

function cluster_synch_error_d_step(x_0, σ, N, eps, a, b, G, uni_clusters; t_final=1000.0) # clusters should be processed
    prob = ODEProblem((dx, x, params, t) -> coupled_fhn_eom!(dx, x, params[1], params[2], params[3], G, b, N), x_0, (0.0, t_final), [a, eps, σ])
    sol = solve(prob; dtmax=0.5)
    cluster_synch_errors = zeros(length(uni_clusters))
    for (i, cluster) in enumerate(uni_clusters)
        t_values, synch_error = local_synch_error(sol, cluster)
        cluster_synch_errors[i] = trapz(t_values, synch_error)/t_final
    end
    global_synch_error = trapz(synch_error_time_series(sol)...)/t_final
    return global_synch_error, cluster_synch_errors, sol.u[end]
end

function cluster_synch_error_sweep(start, stop, init_x0, N_d, N, eps, a, b, G, unique_clusters)
    local d_sweep = range(start, stop, length=N_d)
    local global_synch_error_d_vals = zeros(N_d)
    local cluster_synch_error_d_vals = zeros(N_d, length(unique_clusters))
    x_0 = init_x0
    @showprogress for i in 1:N_d
        global_synch_error_d_vals[i], cluster_synch_error_d_vals[i, :], x_f = cluster_synch_error_d_step(x_0, d_sweep[i], N, eps, a, b, G, unique_clusters)
        x_0 = x_f
    end
    return d_sweep, global_synch_error_d_vals, cluster_synch_error_d_vals
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


N = 10
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)
G = test_matrix_for_cluster_synch() # wattsstrogatzmatrix(N, 3, 1) #

zero_msf = msf_zero()
println("zero_msf = ", zero_msf)

println("msf = ", master_stability_function(zero_msf, 0))
eigenvalues, eigenvectors, clusters, s_matrices = s_matrix_method(G)
unique_clusters = process_clusters(unique(clusters))

eigenvalues = unique(round.(eigenvalues, digits=8))[2:end]
critical_couplings = zero_msf./eigenvalues #unique(zero_msf./eigenvalues)[2:end]#
println("Critical Couplings: ", critical_couplings)
println("Eigenvalues: ", eigenvalues)


N_d = 200
N_realizations = 50
global_forward = zeros(N_realizations, N_d)
global_backward = zeros(N_realizations, N_d)
cluster_forward = zeros(N_realizations, N_d, length(unique_clusters))
cluster_backward = zeros(N_realizations, N_d, length(unique_clusters))

for i in 1:N_realizations
    local G = test_matrix_for_cluster_synch()
    println("Realization ", i)
    global forward_d_sweep, global_forward_d, cluster_forward_d = cluster_synch_error_sweep(0.025, 0.00, zeros(2*N) .+ randn(2*N) .* 0.01, N_d, N, eps, a, b, G, unique_clusters)
    global backward_d_sweep, global_backward_d, cluster_backward_d = cluster_synch_error_sweep(0.00, 0.025, zeros(2*N) .+ randn(2*N) .* 0.01, N_d, N, eps, a, b, G, unique_clusters)
    global_forward[i, :] .= global_forward_d
    global_backward[i, :] .= global_backward_d
    cluster_forward[i, :, :] .= cluster_forward_d
    cluster_backward[i, :, :] .= cluster_backward_d
end

global_forward_avg = mean(global_forward, dims=1)[1, :]
global_backward_avg = mean(global_backward, dims=1)[1, :]
cluster_forward_avg = mean(cluster_forward, dims=1)[1, :, :]
cluster_backward_avg = mean(cluster_backward, dims=1)[1, :, :]

using Serialization
serialize("cluster_sweep_data/global_forward", global_forward)
serialize("cluster_sweep_data/global_backward", global_backward)
serialize("cluster_sweep_data/cluster_forward", cluster_forward)
serialize("cluster_sweep_data/cluster_backward", cluster_backward)
serialize("cluster_sweep_data/forward_d_sweep", forward_d_sweep)
serialize("cluster_sweep_data/backward_d_sweep", backward_d_sweep) # This will do until I write the poster.

f = Figure(size = (700, 900))
ax = Axis(f[1, 1])
ax.xlabel = "Coupling"
ax.ylabel = L"\langle E_{\text{Synch}}\rangle_T"
scatter!(ax, forward_d_sweep[5:end-5], global_forward_avg[5:end-5], label="Right to left") # 5000 timesteps for thermalization
scatter!(ax, backward_d_sweep[5:end-5], global_backward_avg[5:end-5], label="Left to right")
axislegend(position=:rt)
vlines!(ax, critical_couplings[1:3]; label="Critical Couplings", linewidth=1, color = :red)

for i in 1:length(unique_clusters)
    new_ax = Axis(f[i+1, 1])
    new_ax.xlabel = "Coupling"
    new_ax.ylabel = L"\langle E_{\text{Synch}}\rangle_T"
    cluster = unique_clusters[i]
    scatter!(new_ax, forward_d_sweep[6:end-6], cluster_forward_avg[6:end-6, i], label="Cluster $cluster, Right to left") # 6000 timesteps for thermalization
    scatter!(new_ax, backward_d_sweep[6:end-6], cluster_backward_avg[6:end-6, i], label="Cluster $cluster, Left to right")
    vlines!(new_ax, critical_couplings[1:3]; label="Critical Couplings", linewidth=1, color = :red)
    if i == 1
        axislegend(position=:lt)
    else
        axislegend(position=:rt)
    end
end
save("clusters_adiabatic_sweep.png", f)
f

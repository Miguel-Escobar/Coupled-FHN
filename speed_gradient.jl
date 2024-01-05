using DifferentialEquations
using StaticArrays
using Statistics

include("network_simulation.jl")
include("network_matrices_creation.jl")
include("network_sol_analysis.jl")

"""
Calculate the control function for a the coupled neuron system.

Parameters:
    - `gain`: Float64, control gain
    - `each_neuron`: 2D array, state variables of each neuron
    - `coupling_matrix`: 2D array, coupling structure
    - `coupling_terms`: 2D array, coupling terms (given by the coupling Jacobian matrix times the state variables)

Returns:
    - Float64, control signal
"""
function control(gain, each_neuron, coupling_matrix, coupling_terms)
    N = length(coupling_matrix[1, :])
    coupling_strength = 0
    for i in range(1, N)
        for j in range(1, N)
            for k in range(1, N)
                coupling_strength += (each_neuron[1, i] - each_neuron[1, j]) * coupling_matrix[i, k] * coupling_terms[1, k]
                coupling_strength += (each_neuron[2, i] - each_neuron[2, j]) * coupling_matrix[i, k] * coupling_terms[2, k]
            end
        end
    end
    return gain * coupling_strength
end

"""
FitzHugh-Nagumo model equations for a coupled system with control.

Parameters:
    - `dx`: 1D array, output for derivatives
    - `x`: 1D array, state variables
    - `a`: Float64, model parameter
    - `eps`: Float64, model parameter
    - `σ`: Float64, external control signal
    - `control_gain`: Float64, for speed gradient algorithm
    - `coupling_matrix`: 2D array, coupling structure
    - `coupling_jac`: 2D array, coupling Jacobian matrix
"""
function coupled_controlled_fhn_eom!(dx, x, a, eps, σ, control_gain, coupling_matrix, coupling_jac)
    N = length(coupling_matrix[1, :])
    eachneuron = reshape(x, (2, N))
    coupling_terms = coupling_jac * eachneuron
    coupling_strength = control(control_gain, eachneuron, coupling_matrix, coupling_terms) #+ σ
    for i in range(1, N)
        dx_i = fhn_eom(eachneuron[:, i], [a, eps]) .+ coupling_strength .* sum([coupling_matrix[i, j] .* coupling_terms[:, j] for j in 1:N])
        dx[2*i-1:2*i] = dx_i
    end
end


function goal_function(eachneuron)
    goal = 0
    for i in 1:length(eachneuron[1,:])
        for j in 1:length(eachneuron[1, :])
            goal += (eachneuron[1, i] - eachneuron[1, j])^2 + (eachneuron[2, i] - eachneuron[2, j])^2
        end
    end
    return 0.5*goal
end

function control_and_goal_time_series(sol, gain, coupling_matrix, coupling_jac)
    t_values = sol.t
    x_values = sol.u
    control_values = zeros(length(t_values))
    goal_values = zeros(length(t_values))
    for i in 1:length(t_values)
        eachneuron = reshape(x_values[i], (2, N))
        coupling_terms = coupling_jac * eachneuron
        control_values[i] = control(gain, eachneuron, coupling_matrix, coupling_terms)
        goal_values[i] = goal_function(eachneuron)
    end
    return t_values, control_values, goal_values
end


N = 12
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)
σ = 1/N #0.0506 # Coupling strength
γ = abs(σ) # Control gain
G = ring_coupling(N; neighbors=1) # wattsstrogatzmatrix(N, 2, 0.232) #
x_0 = zeros(2*N)
x_0[1:2] .+= 0.1
prob = ODEProblem((dx, x, params, t) -> coupled_fhn_eom!(dx, x, params[1], params[2], params[3], G, b), x_0, (0.0, 25.0), [a, eps, σ])
alg = Tsit5()
sol = solve(prob, alg)

new_x_0 = sol.u[end]
controlled_prob = ODEProblem((dx, x, params, t) -> coupled_controlled_fhn_eom!(dx, x, params[1], params[2], params[3], params[4], G, b), new_x_0, (0.0, 100.0), [a, eps, σ, γ])
controlled_sol = solve(controlled_prob, alg)

uncontrolled_prob = ODEProblem((dx, x, params, t) -> coupled_fhn_eom!(dx, x, params[1], params[2], params[3], G, b), new_x_0, (0.0, 100.0), [a, eps, σ])
uncontrolled_sol = solve(uncontrolled_prob, alg)

using Plots
#plot(new_sol, xlabel="Time", ylabel="System Variables", dpi=600)
controlled_t_val, controlled_kuramoto_val = kuramoto_time_series(controlled_sol, N)
uncontrolled_t_val, uncontrolled_kuramoto_val = kuramoto_time_series(uncontrolled_sol, N)
_, control_values, goal_values = control_and_goal_time_series(controlled_sol, γ, G, b)
l = @layout [a  b ; c]
p1 = plot(uncontrolled_t_val, uncontrolled_kuramoto_val, label="Uncontrolled", xlabel="Time", ylabel="Kuramoto", dpi=600)
plot!(p1, controlled_t_val, controlled_kuramoto_val, label="Controlled", xlabel="Time", ylabel="Kuramoto", dpi=600)
p2 = plot(controlled_t_val, control_values, xlabel="Time", ylabel="Control", dpi=600)
p3 = plot(controlled_t_val, goal_values, xlabel="Time", ylabel="Goal", dpi=600)

observables = plot(p3, p2, p1, layout=l, size=(450, 300), dpi=600)

system_vars = plot(controlled_sol, xlabel="Time", ylabel="System Variables", dpi=600, legend=false)

uncontrolled_vars = plot(uncontrolled_sol, xlabel="Time", ylabel="System Variables", dpi=600, legend=false)

# savefig(observables, "control_observables.png")
# savefig(system_vars, "control_system_vars.png")
# savefig(uncontrolled_vars, "uncontrolled_system_vars.png")


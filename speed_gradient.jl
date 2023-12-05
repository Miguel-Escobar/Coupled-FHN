using DifferentialEquations
using StaticArrays
using Statistics

"""
FitzHugh-Nagumo model equations.

Parameters:
    - `x`: 1D array, state variables [voltage, recovery variable]
    - `params`: 1D array, model parameters [a, epsilon]

Returns:
    - 1D array, derivatives [dv/dt, du/dt]
"""
function fhn_eom(x, params)
    a = params[1]
    eps = params[2]
    dx = (x[1] - (x[1]^3)/3 - x[2])/eps
    dy = x[1] + a
    return [dx, dy]
end

function bmatrix(phi, eps)
    return [cos(phi)/eps sin(phi)/eps; -sin(phi) cos(phi)]
end


"""
FitzHugh-Nagumo model equations for a coupled system without control.

Parameters:
    - `dx`: 1D array, output for derivatives
    - `x`: 1D array, state variables
    - `a`: Float64, model parameter
    - `eps`: Float64, model parameter
    - `coupling_strength`: Float64
    - `coupling_matrix`: 2D array, the connectivity matrix, with the diagonal elements being the negative of the sum of the other elements in the row
    - `coupling_jac`: 2D array, coupling function Jacobian matrix.
"""
function coupled_fhn_eom!(dx, x, a, eps, coupling_strength, coupling_matrix, coupling_jac) # not controlled
    N = length(coupling_matrix[1, :])
    eachneuron = reshape(x, (2, N))
    coupling_terms = coupling_jac * eachneuron
    for i in range(1, N)
        dx_i = fhn_eom(eachneuron[:, i], [a, eps]) .+ coupling_strength .* sum([coupling_matrix[i, j]  .* coupling_terms[:, j] for j in 1:N])
        dx[2*i-1:2*i] = dx_i
    end
end


"""
Calculate the control signal for a coupled system.

Parameters:
    - `gain`: Float64, control gain
    - `each_neuron`: 2D array, state variables of each neuron
    - `coupling_matrix`: 2D array, coupling structure
    - `coupling_terms`: 2D array, coupling terms

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
    return gain * 2 * coupling_strength
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
    coupling_strength = control(control_gain, eachneuron, coupling_matrix, coupling_terms) + σ
    for i in range(1, N)
        dx_i = fhn_eom(eachneuron[:, i], [a, eps]) .+ coupling_strength .* sum([coupling_matrix[i, j]  .* coupling_terms[:, j] for j in 1:N])
        dx[2*i-1:2*i] = dx_i
    end
end

"""
Generate a coupling matrix for a ring network.

Parameters:
    - `size`: Int, number of neurons
    - `neighbors`: Int, number of neighbors to connect to each side (node degree = 2*neighbors)

Returns:
    - 2D array, coupling matrix
"""
function ring_coupling(size; neighbors=1)
    coupling_matrix = zeros(size, size)
    if size > 2*neighbors
        correction = -2*neighbors
    else
        correction = -size + 1
    end
    
    for i in 1:size
        if i + neighbors ≤ size && i - neighbors ≥ 1
            coupling_matrix[i, (i .+ (1:neighbors))] .+= 1
            coupling_matrix[i, i .- (1:neighbors)] .+= 1
            coupling_matrix[i, i] = correction
        else
            indices = unique([mod(j,1:size) for j in (i .- neighbors):(i .+ neighbors) if j != i])
            coupling_matrix[i, indices] .+= 1
            coupling_matrix[i, i] = correction
        end
    end
    return coupling_matrix
end


"""
Generate a Watts-Strogatz small-world network coupling matrix.

Parameters:
    - `size`: Int, number of neurons
    - `neighbors`: Int, number of neighbors to connect initially to each side (node degree = 2*neighbors)
    - `rewiring_prob`: Float64, probability of rewiring an edge

Returns:
    - 2D array, coupling matrix
"""
function wattsstrogatzmatrix(size, neighbors, rewiring_prob)
    coupling_matrix = ring_coupling(size; neighbors=neighbors)
    for i in 1:size
        for j in i:size
            if coupling_matrix[i, j] == 1
                if rand() < rewiring_prob
                    rand_index = rand(1:size)
                    while rand_index == i || coupling_matrix[i, rand_index] == 1
                        rand_index = rand(1:size)
                    end
                    coupling_matrix[i, j] = 0
                    coupling_matrix[i, rand_index] = 1
                    coupling_matrix[rand_index, i] = 1
                end
            end
        end
    end
    return coupling_matrix
end


"""
Calculate the standard deviation of the state variables across neurons.

Parameters:
    - `reshaped_x`: 2D array, reshaped state variables

Returns:
    - Float64, standard deviation
"""
function state_vector_std(reshaped_x)
    return sqrt(var(reshaped_x[1, :]) + var(reshaped_x[2, :]))
end


function std_time_series(sol)
    t_values = sol.t
    x_values = sol.u
    std_values = zeros(length(t_values))
    for i in 1:length(t_values)
        eachneuron = reshape(x_values[i], (2, N))
        std_values[i] = state_vector_std(eachneuron)
    end
    return t_values, std_values
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

function control_and_goal_time_series(sol, coupling_matrix, coupling_jac)
    t_values = sol.t
    x_values = sol.u
    control_values = zeros(length(t_values))
    goal_values = zeros(length(t_values))
    for i in 1:length(t_values)
        eachneuron = reshape(x_values[i], (2, N))
        coupling_terms = coupling_jac * eachneuron
        control_values[i] = control(eachneuron, coupling_matrix, coupling_terms)
        goal_values[i] = goal_function(eachneuron)
    end
    return t_values, control_values, goal_values
end



function kuramoto_time_series(sol, N)
    t_values = sol.t
    x_values = sol.u
    kuramoto_values = zeros(length(t_values))
    for i in 1:length(t_values)
        eachneuron = reshape(x_values[i], (2, N))
        eachangle = atan.(eachneuron[2, :], eachneuron[1, :])
        kuramoto = mean(exp.(im .* eachangle))
        kuramoto_values[i] = abs(kuramoto)
    end
    return t_values, kuramoto_values
end

N = 90
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)
σ = 0.0506 # Coupling strength
γ = 0.01 # Control gain
G = ring_coupling(N; neighbors=2) # wattsstrogatzmatrix(N, 2, 0.232) #
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
_, control_values, goal_values = control_and_goal_time_series(controlled_sol, G, b)
l = @layout [a  b ; c]
p1 = plot(uncontrolled_t_val, uncontrolled_kuramoto_val, label="Uncontrolled", xlabel="Time", ylabel="Kuramoto", dpi=600)
plot!(p1, controlled_t_val, controlled_kuramoto_val, label="Controlled", xlabel="Time", ylabel="Kuramoto", dpi=600)
p2 = plot(controlled_t_val, control_values, xlabel="Time", ylabel="Control", dpi=600)
p3 = plot(controlled_t_val, goal_values, xlabel="Time", ylabel="Goal", dpi=600)

observables = plot(p1, p2, p3, layout=l)

system_vars = plot(controlled_sol, xlabel="Time", ylabel="System Variables", dpi=600, legend=false)

display(observables)
display(system_vars)

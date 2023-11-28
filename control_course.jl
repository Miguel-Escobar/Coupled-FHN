using DifferentialEquations
using StaticArrays
using Statistics

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

function coupled_fhn_eom!(dx, x, a, eps, coupling_strength, coupling_matrix, coupling_jac) # not controlled
    N = length(coupling_matrix[1, :])
    eachneuron = reshape(x, (2, N))
    coupling_terms = coupling_jac * eachneuron
    for i in range(1, N)
        dx_i = fhn_eom(eachneuron[:, i], [a, eps]) .+ coupling_strength .* sum([coupling_matrix[i, j]  .* coupling_terms[:, j] for j in 1:N])
        # for j in 1:N
        #     dx_i .+= coupling_matrix[i, j] .* coupling_strength .* coupling_terms[:, j]
        # end
        dx[2*i-1:2*i] = dx_i
    end
end

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


function wattsstrogatzmatrix(size, neighbors, rewiring_prob)
    coupling_matrix = ring_coupling(size; neighbors=neighbors)
    for i in 1:size
        for j in i:size
            if coupling_matrix[i, j] == 1
                if rand() < rewiring_prob
                    coupling_matrix[i, j] = 0
                    rand_index = rand(1:size)
                    while rand_index == i || coupling_matrix[i, rand_index] == 1
                        rand_index = rand(1:size)
                    end
                    coupling_matrix[i, rand_index] = 1
                    coupling_matrix[rand_index, i] = 1
                end
            end
        end
    end
    return coupling_matrix
end


function state_vector_std(reshaped_x, N)
    return sqrt(var(reshaped_x[1, :]) + var(reshaped_x[2, :]))
end

function std_time_series(sol, N)
    t_values = sol.t
    x_values = sol.u
    std_values = zeros(length(t_values))
    for i in 1:length(t_values)
        eachneuron = reshape(x_values[i], (2, N))
        std_values[i] = state_vector_std(eachneuron, N)
    end
    return t_values, std_values
end

N = 90
eps = 0.05
a = 0.5
b = bmatrix(pi/2-0.1, eps)
σ = 0.0506
G = wattsstrogatzmatrix(N, 3, 0.232)
x_0 = zeros(2*N)
x_0[1:2] .+= 0.1
prob = ODEProblem((dx, x, params, t) -> coupled_fhn_eom!(dx, x, params[1], params[2], params[3], G, b), x_0, (0.0, 100.0), [a, eps, σ])
alg = Tsit5()
sol = solve(prob, alg; dtmax=0.01)
using Plots
# plot(sol, xlabel="Time", ylabel="System Variables", dpi=600)
t_val, std_val = std_time_series(sol, N)
plot(t_val, std_val)
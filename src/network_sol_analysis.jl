using DSP # For the Hilbert transform.
using GLMakie

function state_vector_synch_error(reshaped_x)
    return sqrt(var(reshaped_x[1, :]) + var(reshaped_x[2, :]))
end

function synch_error_time_series(sol)
    t_values = sol.t
    x_values = sol.u
    std_values = zeros(length(t_values))
    for i in 1:length(t_values)
        eachneuron = reshape(x_values[i], (2, N))
        std_values[i] = state_vector_synch_error(eachneuron)
    end
    return t_values, std_values
end

function local_synch_error(sol, neuron_cluster)
    t_values = sol.t
    x_values = sol.u
    local_synch_error_values = zeros(length(t_values))
    for i in 1:length(t_values)
        eachneuron = reshape(x_values[i], (2, N))
        local_synch_error_values[i] = state_vector_synch_error(eachneuron[:, neuron_cluster]) # Test for this.
    end
    return t_values, local_synch_error_values
end


function phase(first_coord)
    phase = angle.(first_coord .+ im .* hilbert(first_coord))
    return phase
end

# function dynamical_phase(t_vals, x_vals, reference_t, reference_x, period)
#     for t in t_vals # Pending

function st_plot(sol)
    t_values = sol.t
    x_values = sol.u
    n_neurons = length(x_values[1, :]) ÷ 2
    u_values = x_values[:, 2 .* (1:n_neurons) .- 1]
    fig = Figure()
    ax = Axis(fig[1, 1])
    heatmap!(ax, t_values, 1:n_neurons, u_values, colormap = :viridis)
    return fig
end

function kuramoto_time_series(sol, N)
    t_values = sol.t
    x_values = hcat(sol.u...)
    kuramoto_values = zeros(length(t_values))
    for i in 1:N
        first_coord = x_values[2*i-1, :]
        kuramoto_values += exp.(im .* phase(first_coord))
    end
    kuramoto_values = abs.(kuramoto_values ./ N)
    return t_values, kuramoto_values
end

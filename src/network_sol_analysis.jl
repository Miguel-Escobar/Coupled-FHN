using DSP

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


function phase(first_coord)
    phase = angle.(first_coord .+ im .* hilbert(first_coord))
    return phase
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
    # kuramoto_values = zeros(length(t_values))
    # for i in 1:length(t_values)
    #     eachneuron = reshape(x_values[i], (2, N))
    #     eachangle = atan.(eachneuron[2, :], eachneuron[1, :])
    #     kuramoto = mean(exp.(im .* eachangle))
    #     kuramoto_values[i] = abs(kuramoto)
    # end
    # return t_values, kuramoto_values
end

# function local_synch_error(sol, neuron_cluster)

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
                    rand_index = rand(1:size)
                    while rand_index == i || coupling_matrix[i, rand_index] == 1
                        rand_index = rand(1:size)
                    end
                    coupling_matrix[i, j] = 0
                    coupling_matrix[j, i] = 0
                    coupling_matrix[i, rand_index] = 1
                    coupling_matrix[rand_index, i] = 1
                end
            end
        end
    end
    return coupling_matrix
end
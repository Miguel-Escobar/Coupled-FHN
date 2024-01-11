using Graphs
using LinearAlgebra
using StaticArrays

function test_matrix_for_cluster_synch()
    [9 -1 -1 -1 -1 -1 -1 -1 -1 -1;
    -1 9 -1 -1 -1 -1 -1 -1 -1 -1;
    -1 -1 9 -1 -1 -1 -1 -1 -1 -1;
    -1 -1 -1 15 -2 -2 -2 -2 -2 -2;
    -1 -1 -1 -2 15 -2 -2 -2 -2 -2;
    -1 -1 -1 -2 -2 15 -2 -2 -2 -2;
    -1 -1 -1 -2 -2 -2 18 -3 -3 -3;
    -1 -1 -1 -2 -2 -2 -3 18 -3 -3;
    -1 -1 -1 -2 -2 -2 -3 -3 18 -3;
    -1 -1 -1 -2 -2 -2 -3 -3 -3 18]
end

function complete_network(size)
    coupling_matrix = zeros(size, size)
    coupling_matrix .= 1
    coupling_matrix = coupling_matrix - Diagonal(vec(sum(coupling_matrix, dims=2))) # Ensure zero row sum
    return -coupling_matrix
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
    return -coupling_matrix
end

function wattsstrogatzmatrix(size, neighbors, rewiring_prob)
    g = watts_strogatz(size, 2*neighbors, rewiring_prob)
    coupling_matrix = adjacency_matrix(g) # Transform into full matrix (not sparse) (will test)
    coupling_matrix = Matrix(coupling_matrix)
    coupling_matrix = coupling_matrix - Diagonal(vec(sum(coupling_matrix, dims=2))) # Ensure zero row sum
    return -coupling_matrix
end
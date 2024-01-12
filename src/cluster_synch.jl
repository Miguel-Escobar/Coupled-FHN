# Following a recipe by A. Bayani et. al.
using LinearAlgebra
using Graphs

function create_S_matrices(eigenvectors)
    returnable = []
    n_eigvectors = length(eigenvectors[1, :])
    S = zeros(n_eigvectors, n_eigvectors)
    for i in reverse(1:n_eigvectors)
        E = zeros(n_eigvectors, n_eigvectors)
        for j in 1:n_eigvectors
            for k in 1:n_eigvectors
                E[j, k] = (eigenvectors[j, i] - eigenvectors[k, i])^2
            end
        end
        S += E
        push!(returnable, S)
    end
    return reverse(returnable)
end

function s_matrix_method(matrix)
    eigenvalues, eigenvectors = eigen(matrix)
    s_matrices = create_S_matrices(eigenvectors)
    clusters = []
    for s_matrix in s_matrices
        cluster_indices = isapprox.(s_matrix, 2; atol=1e-8)
        graph = SimpleGraph(cluster_indices)
        connected_comp = [comp for comp in connected_components(graph) if length(comp) > 1]
        push!(clusters, connected_comp)
    end
    return eigenvalues, eigenvectors, clusters, s_matrices
end
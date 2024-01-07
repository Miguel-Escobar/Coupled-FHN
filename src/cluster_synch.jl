# Following a recipe by A. Bayani et. al.
using LinearAlgebra

function create_S_matrices(eigenvectors)
    returnable = []
    n_eigvectors = length(eigenvectors[1, :])
    S = zeros(n_eigvectors, n_eigvectors)
    for i in 1:n_eigvectors
        E = zeros(n_eigvectors, n_eigvectors)
        for j in 1:n_eigvectors
            for k in 1:n_eigvectors
                E[j, k] = (eigenvectors[j, i] - eigenvectors[k, i])^2
            end
        end
        S += E
        push!(returnable, S)
    end
    return returnable
end

function s_matrix_method(matrix)
    eigenvalues, eigenvectors = eigen(matrix)
    s_matrices = create_S_matrices(eigenvectors)
    clusters = [] 
    for i in 1:length(s_matrices)
        cluster_indices = findall(x -> isapprox(x, 2; atol=1e-8), s_matrices[i])
        push!(clusters, cluster_indices)
    end
    return eigenvalues, eigenvectors, clusters, s_matrices
end

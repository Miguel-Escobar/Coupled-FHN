using DifferentialEquations # For later, so that I don't have to import it when I use this code.
using StaticArrays
using Tullio: @tullio
using LoopVectorization

function fhn_eom(x, params)
    a = params[1]
    eps = params[2]
    dx = (x[1] - (x[1]^3)/3 - x[2])/eps
    dy = x[1] + a
    return SVector{2}(dx, dy)
end

function bmatrix(phi, eps)
    return -[cos(phi)/eps sin(phi)/eps; -sin(phi) cos(phi)]
end

function coupled_fhn_eom!(dx::Vector{Float64}, x::Vector{Float64}, a::Float64, eps::Float64, coupling_strength::Float64, coupling_matrix::Matrix{Float64}, b::Matrix{Float64}, N::Int64)
    eachneuron = reshape(x, (2, N))
    coupling_terms = b * eachneuron
    coupling = zeros(Float64, N, 2)
    @tullio coupling[i, j] = coupling_matrix[i, r] * coupling_terms[j, r]
    eom_params = SVector{2}(a, eps)
    for i in range(1, N)
        thisneuron = @view eachneuron[:, i]
        thiscoupling = @view coupling[i, :]
        dx[2*i-1:2*i] .= fhn_eom(thisneuron, eom_params) .+ coupling_strength .* thiscoupling
    end
    nothing
end

coupled_fhn_eom!(dx::Vector{Float64}, x::Vector{Float64}, a::Float64, eps::Float64, coupling_strength::Float64, coupling_matrix::Matrix{Int64}, b::Matrix{Float64}, N::Int64) = coupled_fhn_eom!(dx, x, a, eps, coupling_strength, convert(Matrix{Float64}, coupling_matrix), b, N)
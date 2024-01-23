using StaticArrays, LinearAlgebra
using DynamicalSystems
using OrdinaryDiffEq
using LaTeXStrings
using GLMakie
using ProgressMeter
using Base.Threads
using Roots

function fhn_eom(x, params, t)
    a = params[1]
    eps = params[2]
    dx = (x[1] - (x[1]^3)/3 - x[2])/eps
    dy = x[1] + a
    return SVector(dx, dy)
end

function fhn_jac(x, params, t)
    eps = params[2]
    dx_dx = (1 - x[1]^2)/eps
    dx_dy = -1/eps
    dy_dx = 1
    dy_dy = 0
    returnable = SMatrix{2, 2, Float64}(dx_dx, dy_dx, dx_dy, dy_dy)  #SA_F64[dx_dx dx_dy; dy_dx dy_dy] 
    return returnable
end

function msf_eom!(dxchi, xchi, params, t)
    aeps = @view params[1:2]
    alpha = params[3]
    beta = params[4]
    c = params[5] 
    phi = params[6]
    B = couplingJacobian(phi, aeps[2])
    chireal = @view xchi[3:4]
    chiimag = @view xchi[5:6]
    x = @view xchi[1:2]
    jac = fhn_jac(x, aeps, t)
    dxchi[1:2] .= fhn_eom(x, aeps, t)
    dxchi[3:4] .= jac * chireal .+ c.*(alpha.*B*chireal .- beta.*B*chiimag)
    dxchi[5:6] .= jac * chiimag .+ c.*(alpha.*B*chiimag .+ beta.*B*chireal)
    return nothing
end

function couplingJacobian(phi, eps)
    #return -1 .* SA[cos(phi)/eps sin(phi)/eps; -sin(phi) cos(phi)]
    return -1 .* SMatrix{2, 2, Float64}(cos(phi)/eps, -sin(phi), sin(phi)/eps, cos(phi))
end

function msf_system(alpha, beta; a=0.5, eps=0.05, coupling=1.0, phi=(pi/2)-0.1, diffeq=(alg=Tsit5(), abstol = 1e-9, reltol = 1e-9))
    ds = ContinuousDynamicalSystem(msf_eom!, SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], SA[a, eps, alpha, beta, coupling, phi], diffeq=diffeq)
    return ds
end

function master_stability_function(alpha, beta; testfunc=(state1, d0) -> [state1[1:2] ; state1[3:end] .+ d0/sqrt(4)], kwargs...)
    system = msf_system(alpha, beta; kwargs...)
    return lyapunov(system, 1000.0; Δt = 0.1, Ttr=100.0, inittest=testfunc, d0=1e-9)
end

function plot_msf_regions(n_rows; kwargs...)
    alpha_sweep = range(-1.5, 1.5, length=n_rows)
    beta_sweep = range(-0.5, 0.5, length=n_rows)
    msf = zeros(length(alpha_sweep), length(beta_sweep))

    @showprogress for i in 1:length(alpha_sweep)
        alpha = alpha_sweep[i]
        Threads.@threads for j in 1:length(beta_sweep)
            beta = beta_sweep[j]
            msf[i, j] = master_stability_function(alpha, beta; kwargs...)
        end
    end

    levels = [-1e10, 0, 1e10]
    println(alpha_sweep)
    println(beta_sweep)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig, xlabel=L"α", ylabel=L"β", zlabel=L"λ")
    contour!(ax, alpha_sweep, beta_sweep, msf, levels=levels, fill=true)
    display(fig)
end

function plot_msf_vs_eigs(start, stop, n_points; kwargs...)
    eigenvalue_real_sweep = range(start, stop, length=n_points)
    msf_sweep = zeros(length(eigenvalue_real_sweep))

    Threads.@threads for i in 1:length(eigenvalue_real_sweep)
        msf_sweep[i] = master_stability_function(eigenvalue_real_sweep[i], 0.0; kwargs...)
    end
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig)
    lines!(ax, eigenvalue_real_sweep, msf_sweep)
    display(fig)
end

function synch_is_stable(coupling_matrix; kwargs...)
    eigenvalues = eigvals(coupling_matrix)
    msf_for_eigs = master_stability_function.(real.(eigenvalues), imag.(eigenvalues); kwargs...)
    return all(msf_for_eigs .≤ 0)
end


function plot_msf_regions_with_eigs(n_rows, coupling_matrix; savefigure=false, kwargs...)
    alpha_sweep = range(-1.5, 0.5, length=n_rows)
    beta_sweep = range(-0.5, 0.5, length=n_rows)
    msf = zeros(length(alpha_sweep), length(beta_sweep))

    @showprogress for j in 1:length(alpha_sweep)
        alpha = alpha_sweep[j]
        Threads.@threads for i in 1:length(beta_sweep)
            beta = beta_sweep[i]
            msf[i, j] = master_stability_function(alpha, beta; kwargs...)
        end
    end

    eigs = eigvals(coupling_matrix)
    msf_for_eigs = master_stability_function.(real.(eigs), imag.(eigs); kwargs...)
    levels = [-1e10, 0, 1e10]

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig, xlabel=L"α", ylabel=L"β", zlabel=L"λ", xlims=(-1.5, 0.5), ylims=(-0.5, 0.5))
    contour!(ax, alpha_sweep, beta_sweep, msf, levels=levels, fill=true)
    scatter!(ax, real.(eigs), imag.(eigs), color=:red, label="Eigenvalues")
    display(fig)

    if savefigure
        save("msf_with_eigs.png", fig)
    end

    if all(msf_for_eigs .≤ 0)
        println("Synchronization is stable")
        println("Max Lyapunov: ", maximum(msf_for_eigs))
    else
        println("Synchronization is unstable")
        println("Max Lyapunov: ", maximum(msf_for_eigs))
    end
end

function msf_zero(kwargs...)
    return find_zero(alpha -> master_stability_function(alpha, 0.0; kwargs...), 0.2)
end

using Pkg

dependencies = ["DifferentialEquations",
                "OrdinaryDiffEq", # For solver options in msf.jl I think.
                "StaticArrays",
                "DynamicalSystems",
                "Plots",
                "LaTeXStrings",
                "ProgressMeter",
                "Peaks",
                "Graphs",
                "GLMakie"]

for dep in dependencies
    Pkg.add(dep)
end

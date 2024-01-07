using Pkg

dependencies = ["DifferentialEquations",
                "OrdinaryDiffEq",
                "StaticArrays",
                "DynamicalSystems",
                "Plots",
                "LaTeXStrings",
                "ProgressMeter",
                "Peaks",
                "Graphs"]

for dep in dependencies
    Pkg.add(dep)
end

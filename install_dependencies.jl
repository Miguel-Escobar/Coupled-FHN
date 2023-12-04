using Pkg

dependencies = ["DifferentialEquations",
                "OrdinaryDiffEq",
                "StaticArrays",
                "DynamicalSystems",
                "Plots",
                "LaTeXStrings",
                "ProgressMeter",
                "Peaks"]

for dep in dependencies
    Pkg.add(dep)
end

module SMeshSmokeValidation

# This file is formatted to match StanfordAA228V.jl/src/StanfordAA228V.jl

# External dependencies
using LinearAlgebra
using LazySets
using Plots

# Internal dependencies
# include("Counted.jl")
# using .Counted

# Exported functions
export 
    burn_scene,
    make_scene_from_csv,
    plot_scene,
    animate_burn_scene

# Included files
include("make_scene.jl")
include("plot_scene.jl")

end

module SMeshSmokeValidation

# This file is formatted to match StanfordAA228V.jl/src/StanfordAA228V.jl

# External dependencies
using LinearAlgebra
using LazySets
using Plots
using Images

# Internal dependencies
# include("Counted.jl")
# using .Counted

# Exported functions
export 
    burn_scene,
    parse_coord_str_to_poly,
    make_scene_from_csv,
    sample_smoke_from_poly,
    sample_smoke_from_scene,
    plot_scene,
    animate_burn_scene

# Included files
include("make_scene.jl")
include("plot_scene.jl")
include("smoke.jl")

end

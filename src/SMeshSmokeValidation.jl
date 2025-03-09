module SMeshSmokeValidation

# This file is formatted to match StanfordAA228V.jl/src/StanfordAA228V.jl

# External dependencies
println("Loading External Dependencies...")
using LinearAlgebra
using LazySets
using Plots
using Images

# Internal dependencies
# include("Counted.jl")
# using .Counted

println("Loading SMeshSmokeValidation...")

# Exported functions
export 
    # Make Scene
    BurnScene,
    parse_coord_str_to_poly,
    make_scene_from_csv,
    # Plot Scene
    plot_scene,
    animate_burn_scene,
    # Plume Model
    PlumeFromPointSource,
    query_plume_model,
    # Smoke Points
    sample_smoke_from_poly,
    sample_smoke_from_scene

# Included files
include("make_scene.jl")
include("plot_scene.jl")
include("plume_model.jl")
include("smoke_points.jl")

end

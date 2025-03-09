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

println("Setting exports")

# Exported functions
export 
    # Make Scene
    BurnScene,
    parse_coord_str_to_poly,
    make_scene_from_csv,
    # Plot Scene
    plot_scene,
    animate_burn_scene,
    plot_single_plume_bounds,
    plot_multiple_plumes_bounds,
    # Plume Model
    PlumeFromPointSource,
    query_plume_model,
    # Smoke Points
    sample_smoke_from_poly,
    sample_smoke_from_scene

println("Loading SMeshSmokeValidation...")
# Included files (order needs to match dependencies)
include("make_scene.jl")
include("plume_model.jl")
include("smoke_points.jl")
include("plot_scene.jl")

end

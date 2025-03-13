module SMeshSmokeValidation

# This file is formatted to match StanfordAA228V.jl/src/StanfordAA228V.jl

# External dependencies
println("Loading External Dependencies...")
using LinearAlgebra
using LazySets
using Plots
using Images
using Random
using Distributions

# Internal dependencies
# include("Counted.jl")
# using .Counted

println("Setting Exports")

# Exported functions
export 
    # Height Distribution
    HeightDistribution,
    # Wind Distribution
    WindDistribution,
    # Disturbances
    FullDisturbances,
    FullDisturbanceTrajectory,
    sample_disturbances,
    disturbance_trajectory_likelihood,
    disturbance_trajectory_log_likelihood,
    most_likely_disturbance_trajectories,
    # Make Scene
    BurnScene,
    parse_coord_str_to_poly,
    # make_scene_from_csv, # Depracated
    # load_burn_area_file, # Depracated
    load_burn_scene_from_files,
    # Plume Model
    PlumeFromPointSource,
    query_plume_model,
    query_multiple_plumes_points,
    # Smoke Points
    sample_smoke_from_poly,
    sample_smoke_from_scene,
    # Run Scene
    gen_smoke_samples_over_time,
    gen_plumes_over_time,
    gen_sensor_readings_of_plume,
    # Plot Scene
    plot_scene,
    animate_burn_scene,
    plot_single_plume_bounds,
    plot_multiple_plumes_bounds

println("Loading SMeshSmokeValidation...")
# Included files (order needs to match dependencies)
include("height_distribution.jl")
include("wind_distribution.jl")
include("disturbance.jl")
include("make_scene.jl")
include("plume_model.jl")
include("smoke_points.jl")
include("run_scene.jl")
include("plot_scene.jl")

end

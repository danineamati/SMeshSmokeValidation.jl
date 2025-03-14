###############################
# Load files and import libraries
###############################
using Plots
using Random
using LazySets         # For operations with geometric sets (used in burn scenes)
using Images           # For image loading and processing
using FileIO           # For file I/O related to images
using Colors           # For handling colors
using CSV, DataFrames  # For reading and processing CSV data

using SMeshSmokeValidation

# Set a global RNG seed for reproducibility
Random.seed!(42)

###############################
# Simulation Helper Functions 
###############################

function initialize_scene(dataset::String, simulation::Int, moisture_level::String)
    # Define file paths using joinpath
    burn_area_file       = joinpath("data", "BurnData", dataset, "BurnAreas", "$(dataset)BurnAreas.txt")
    snode_locations_file = joinpath("data", "BurnData", dataset, "SNodeLocations", "snodelocations_sim$(simulation).txt")
    reference_bounds_file= joinpath("data", "BurnData", dataset, "ReferencePoints", "referencelocations.txt")
    total_burn_area_file = joinpath("data", "BurnData", dataset, "BurnAreas", "TotalBurnArea.txt")
    acreage_consumption_file  = joinpath("data", "FuelConsumptionProfiles", "$(dataset)_Scenarios.csv")
    
    # Load the BurnScene using the helper from make_scene.jl.
    burn_scene = load_burn_scene_from_files(burn_area_file, snode_locations_file, reference_bounds_file, total_burn_area_file, acreage_consumption_file, moisture_level)

    # Print details to verify instantiation
    println("Initialized BurnScene:")
    println(" - Burn polygons: ", length(burn_scene.burn_polys_t))
    println(" - Time steps: ", length(burn_scene.t))
    println(" - Sensor nodes: ", length(burn_scene.snode_locations))
    println(" - X bounds: ", burn_scene.reference_bounds_x)
    println(" - Y bounds: ", burn_scene.reference_bounds_y)
    println(" - Perimeter sample points: ", length(burn_scene.perimeter_sample_points))
    println(" - Q value: ", burn_scene.Q)
    
    return burn_scene
end

function run_disturbance_example(burn_scene::BurnScene)
    println("\n=== Running Disturbance Example ===")
    
    # Define the save directory using joinpath
    save_dir = joinpath("plots", "disturbance_examples", "Malibu_fuzzed")
    if !isdir(save_dir)
        mkpath(save_dir)
    end

    Q = burn_scene.Q
    # Define disturbance parameter distributions
    default_disturb = FullDisturbances(
        HeightDistribution(1000.0, 80.0),
        WindDistribution(
            [0.0, -120, -180],
            [9.0, 5.0, 5.0],
            [0.7, 0.15, 0.15],
            5.0, 1.0;
            in_deg=true
        )
    )
    malibu_nominal_disturb = FullDisturbances(
        HeightDistribution(Q, 30.0),
        WindDistribution(
            [-120.0, 90.0, 150.0],
            [9.0, 5.0, 5.0],
            [0.1, 0.7, 0.2],
            10.0, 5.0;
            in_deg=true
        )
    )
    malibu_fuzzed_disturb = FullDisturbances(
        HeightDistribution(Q*2, 5.0),
        WindDistribution(
            [-120.0, 90.0, 150.0],
            [9.0, 5.0, 5.0],
            [0.7, 0.15, 0.15],
            5.0, 1.0;
            in_deg=true
        )
    )
    
    # Choose disturbance distributions based on save_dir
    if save_dir == joinpath("plots", "disturbance_examples", "full_disturbance_testing")
        disturb = default_disturb
        disturb_nominal = default_disturb
    elseif save_dir == joinpath("plots", "disturbance_examples", "Malibu_nominal")
        disturb = malibu_nominal_disturb
        disturb_nominal = malibu_nominal_disturb
    elseif save_dir == joinpath("plots", "disturbance_examples", "Malibu_fuzzed")
        disturb = malibu_fuzzed_disturb
        disturb_nominal = malibu_nominal_disturb
    end

    # Sample trajectories and compute likelihoods
    num_timesteps = 10
    num_trajectories = 1000
    trajectories = [ sample_disturbances(disturb, num_timesteps) for _ in 1:num_trajectories ]
    likelihoods_sampled = [ disturbance_trajectory_likelihood(disturb, traj) for traj in trajectories ]
    likelihoods_nominal = [ disturbance_trajectory_likelihood(disturb_nominal, traj) for traj in trajectories ]
    loglikelihoods_sampled = [ disturbance_trajectory_log_likelihood(disturb, traj) for traj in trajectories ]
    loglikelihoods_nominal = [ disturbance_trajectory_log_likelihood(disturb_nominal, traj) for traj in trajectories ]

    # Plot the log-likelihood histograms
    lowest_loglikelihood = min(minimum(loglikelihoods_sampled), minimum(loglikelihoods_nominal))
    highest_loglikelihood = max(maximum(loglikelihoods_sampled), maximum(loglikelihoods_nominal))
    n_bins = 200
    bin_edges = range(lowest_loglikelihood, stop=highest_loglikelihood, length=n_bins+1)
    sampled_label = occursin("fuzzed", lowercase(save_dir)) ? "Fuzzed" : "Sampled"

    h = histogram(loglikelihoods_nominal, label="Nominal",
        alpha=0.5, bins=bin_edges, normalize=:pdf,
        xlabel="Disturbance Trajectory Log Likelihood", ylabel="Estimated PDF",
        framestyle=:box, dpi=300)
    histogram!(loglikelihoods_sampled, label=sampled_label,
        alpha=0.5, bins=bin_edges, normalize=:pdf)
    savefig(h, joinpath(save_dir, "disturbance_loglikelihood_histogram.png"))

    # Identify and display the top k trajectories
    k_top = 10
    top_trajectories, top_indices = most_likely_disturbance_trajectories(loglikelihoods_nominal, trajectories, k_top)
    println("\n--- Top $k_top Disturbance Trajectories (Nominal) ---")
    for i in 1:k_top
        println("Trajectory $i: Log-likelihood = $(loglikelihoods_nominal[top_indices[i]])")
        println("Trajectory $i: ", top_trajectories[i])
    end

    println("=== Disturbance Example Completed ===\n")
end


# Instead of reloading the scene data, we now pass the BurnScene object.
function run_burn_scene_with_background(burn_scene::BurnScene, dataset::String, simulation::Int, moisture_level::String)
    println("\n=== Running Burn Scene with Background ===")
    
    # Set up save directory
    save_dir = joinpath("plots", dataset)
    if !isdir(save_dir)
        mkpath(save_dir)
    end

    # Load the background image
    img_path = joinpath("data", "BurnData", "DEMs_and_Buffered_Burns", "Background_$(dataset).png")
    background_img = load(img_path)
    background_img = RGB.(background_img)
    println("Background image dimensions: ", size(background_img))

    # Create coordinate arrays for the background image
    xlims = burn_scene.reference_bounds_x
    ylims = burn_scene.reference_bounds_y
    nx = size(background_img, 2)
    ny = size(background_img, 1)
    background_x = range(xlims[1], stop=xlims[2], length=nx)
    background_y = range(ylims[1], stop=ylims[2], length=ny)
    
    # Plot the burn scene over time
    for t_ind in 1:(length(burn_scene.t) + 1)
        p = plot_scene(burn_scene, t_ind,
                       n_smoke_samples=10,
                       background_image=background_img,
                       background_x=background_x,
                       background_y=background_y)
        savefig(p, joinpath(save_dir, "burn_scene_test_$(dataset)_topobkgd_$(t_ind).png"))
    end
    
    println("=== Burn Scene with Background Completed ===\n")
end

function run_changing_wind_scene(burn_scene::BurnScene, dataset::String, moisture_level::String)
    println("\n=== Running Changing Wind Scene ===")
    
    # Set up save directory
    save_dir = joinpath("plots", dataset, "changing_wind")
    if !isdir(save_dir)
        mkpath(save_dir)
    end
    
    # Generate smoke samples and plumes using the preloaded scene
    n_smoke_samples = 10
    smoke_samples_per_time = gen_smoke_samples_over_time(burn_scene, n_smoke_samples=n_smoke_samples)
    
    # Use the preloaded Q (or compute it if not present)
    Q = burn_scene.Q
    
    h_plume = 10.0
    air_class = "D"
    plumes_per_time = gen_plumes_over_time(smoke_samples_per_time, Q, h_plume, air_class)
    
    # Use the perimeter pseudo-nodes as sensor locations
    sensor_locations = burn_scene.snode_locations
    sensor_locations = [ [pt[1], pt[2], 0.0] for pt in sensor_locations ]  # extend to 3D if needed
    
    # Define wind vectors
    wind_speed = 10.0
    wind_dirs = [30.0, 35.0, 40.0, 35.0, -30.0, -25.0, -20.0, 35.0, 30.0]
    to_wind_vec(wdir) = wind_speed .* [cosd(wdir), sind(wdir), 0.0]
    wind_vecs = to_wind_vec.(wind_dirs)
    wind_len_mult = 3.0
    vmin = -10.0
    vmax = 0.0
    
    # Compute sensor readings and take log10
    sensor_readings_per_time = gen_sensor_readings_of_plume(plumes_per_time, sensor_locations, wind_vecs)
    sensor_readings_per_time_log = [ log10.(sr) for sr in sensor_readings_per_time ]

    # Compute perimeter readings using perimeter pseudo‑nodes
    perimeter_nodes = burn_scene.perimeter_sample_points
    perimeter_nodes_3D = [ [pt[1], pt[2], 0.0] for pt in perimeter_nodes ]
    perimeter_readings_per_time = gen_sensor_readings_of_plume(plumes_per_time, perimeter_nodes_3D, wind_vecs)
    perimeter_readings_per_time_log = [ log10.(sr) for sr in perimeter_readings_per_time ]

    # Save directories (adjust as needed)
    sensor_save_dir = joinpath(save_dir, "sensor_readings")
    perimeter_save_dir = joinpath(save_dir, "perimeter_readings")
    if !isdir(sensor_save_dir)
        mkpath(sensor_save_dir)
    end
    if !isdir(perimeter_save_dir)
        mkpath(perimeter_save_dir)
    end

    # Annotate sensor readings
    vert_offset = 30
    txt_size = 8

    # --- Plot Sensor Readings ---
    for t in 1:(length(burn_scene.t) + 1)
        p = plot(framestyle=:box)
        curr_wind = wind_vecs[t]
        
        # Plot burn polygons (same for both sets)
        for ind in eachindex(burn_scene.burn_polys_t)
            if ind < t
                plot!(burn_scene.burn_polys_t[ind], color=:black, alpha=0.5, label="", linealpha=0.0)
            elseif ind == t
                plot!(burn_scene.burn_polys_t[ind], color=:red, label="")
            else
                plot!(burn_scene.burn_polys_t[ind], color=:lightgray, alpha=0.5, label="", linealpha=0.0)
            end
        end
        
        # Plot original sensor locations with readings
        scatter!([pt[1] for pt in sensor_locations],
                [pt[2] for pt in sensor_locations],
                zcolor=[sr for sr in sensor_readings_per_time_log[t]],
                label="",
                c=cgrad(:bilbao, rev=true),
                clim=(vmin, vmax))
        

        for (i, pt) in enumerate(sensor_locations)
            annotate!(pt[1], pt[2] + vert_offset,
                    text(round(sensor_readings_per_time_log[t][i], digits=1), txt_size),
                    halign=:center, valign=:center, color="black")
        end
        
        xlims!(burn_scene.reference_bounds_x...)
        ylims!(burn_scene.reference_bounds_y...)
        curr_wind_dir = Int(wind_dirs[t])
        savefig(p, joinpath(sensor_save_dir, "example_sensor_readings_$(t)_$(curr_wind_dir).png"))
    end

    # --- Plot Perimeter Readings ---
    for t in 1:(length(burn_scene.t) + 1)
        p = plot(framestyle=:box)
        curr_wind = wind_vecs[t]
        
        # Plot burn polygons (same as before)
        for ind in eachindex(burn_scene.burn_polys_t)
            if ind < t
                plot!(burn_scene.burn_polys_t[ind], color=:black, alpha=0.5, label="", linealpha=0.0)
            elseif ind == t
                plot!(burn_scene.burn_polys_t[ind], color=:red, label="")
            else
                plot!(burn_scene.burn_polys_t[ind], color=:lightgray, alpha=0.5, label="", linealpha=0.0)
            end
        end
        
        # Plot perimeter sensor locations with readings
        scatter!([pt[1] for pt in perimeter_nodes_3D],
                [pt[2] for pt in perimeter_nodes_3D],
                zcolor=[sr for sr in perimeter_readings_per_time_log[t]],
                label="",
                c=cgrad(:bilbao, rev=true),
                clim=(vmin, vmax))
        
        # Annotate perimeter readings
        for (i, pt) in enumerate(perimeter_nodes_3D)
            annotate!(pt[1], pt[2] + vert_offset,
                    text(round(perimeter_readings_per_time_log[t][i], digits=1), txt_size),
                    halign=:center, valign=:center, color="black")
        end
        
        xlims!(burn_scene.reference_bounds_x...)
        ylims!(burn_scene.reference_bounds_y...)
        curr_wind_dir = Int(wind_dirs[t])
        savefig(p, joinpath(perimeter_save_dir, "example_perimeter_readings_$(t)_$(curr_wind_dir).png"))
    end

    
    # Plot plume maps over time
    println("Plotting plume maps...")
    bounds = vcat(burn_scene.reference_bounds_x..., burn_scene.reference_bounds_y...)
    for t in 1:(length(burn_scene.t) + 1)
        curr_wind = wind_vecs[t]
        p = plot_multiple_plumes_bounds(plumes_per_time[t], bounds, curr_wind,
                                        x_num=101, y_num=103, vmax=vmax)
        if length(smoke_samples_per_time[t]) > 0
            for source_location in eachcol(smoke_samples_per_time[t])
                source_location = source_location[:]  # convert to 1D array
                scatter!([source_location[1]], [source_location[2]], color="gray", label="")
                quiver!([source_location[1]], [source_location[2]],
                        quiver=([wind_len_mult * curr_wind[1]], [wind_len_mult * curr_wind[2]]),
                        color="black", lw=2, label="")
            end
        end
        xlims!(burn_scene.reference_bounds_x...)
        ylims!(burn_scene.reference_bounds_y...)
        curr_wind_dir = Int(wind_dirs[t])
        savefig(p, joinpath(save_dir, "example_plume_multiple_$(t)_$(curr_wind_dir).png"))
    end
    
    println("=== Changing Wind Scene Completed ===\n")
end

###############################
# Main Function
###############################
function main()
    # Define parameters (using the same dataset for both scenes)
    dataset = "Shasta"
    simulation = 1
    moisture_level = "moist"
    
    # Initialize the scene once
    burn_scene = initialize_scene(dataset, simulation, moisture_level)

    run_disturbance_example(burn_scene)
    run_burn_scene_with_background(burn_scene, dataset, simulation, moisture_level)
    run_changing_wind_scene(burn_scene, dataset, moisture_level)
end

main()

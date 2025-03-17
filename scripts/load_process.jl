#######################################
# Load files and import libraries
#######################################
using Plots
using Random
using LazySets         # For operations with geometric sets (used in burn scenes)
using Images           # For image loading and processing
using FileIO           # For file I/O related to images
using Colors           # For handling colors
using CSV, DataFrames  # For reading and processing CSV data
using Distributions    # For sampling from distributions

using SMeshSmokeValidation

# Set a global RNG seed for reproducibility
Random.seed!(42)


# Pull in the wind data
include("wind_configurations.jl")


#######################################
# Check failure 
#######################################

function is_failure(sensor_readings::Vector{Float64}, 
                    perimeter_readings::Vector{Float64};
                    convert_to_linear::Bool=false)
    # Filter all the infinite values
    sensor_readings = filter(x -> isfinite(x), sensor_readings)
    perimeter_readings = filter(x -> isfinite(x), perimeter_readings)

    if convert_to_linear
        sensor_readings = 10 .^ sensor_readings
        perimeter_readings = 10 .^ perimeter_readings
    end

    # If the mean of the perimeter is greater than the mean of the sensor 
    # readings, then the sensor readings are failing to detect the plume
    return Distributions.mean(perimeter_readings) > Distributions.mean(sensor_readings)
    
end



#######################################
# Simulation Helper Functions 
#######################################

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

function run_disturbance_example(burn_scene::BurnScene, dataset::String)
    println("\n=== Running Disturbance Example ===")
    
    # Define the save directory using joinpath
    save_dir = joinpath("plots", "disturbance_examples", dataset)
    if !isdir(save_dir)
        mkpath(save_dir)
    end

    Q = burn_scene.Q
    # Define disturbance parameter distributions
    
    # Choose disturbance distributions based on save_dir
    if dataset == "HenryCoe"
        disturb = HENRY_COE_FUZZED(Q)
        disturb_nominal = HENRY_COE_NOMINAL(Q)
    elseif dataset == "Malibu"
        disturb = MALIBU_FUZZED(Q)
        disturb_nominal = MALIBU_NOMINAL(Q)
    elseif dataset == "Shasta"
        disturb = SHASTA_FUZZED(Q)
        disturb_nominal = SHASTA_NOMINAL(Q)
    else
        println("Dataset not recognized. Using default disturbance.")
        disturb = DEFAULT_DISTURBANCE(Q)
        disturb_nominal = DEFAULT_DISTURBANCE(Q)
    end

    # Sample trajectories and compute likelihoods
    num_timesteps = 10
    num_trajectories = 1000
    trajectories = [ sample_disturbances(disturb, num_timesteps) for _ in 1:num_trajectories ]
    # We will use the log-likelihoods for numerical stability
    # likelihoods_sampled = [ disturbance_trajectory_likelihood(disturb, traj) for traj in trajectories ]
    # likelihoods_nominal = [ disturbance_trajectory_likelihood(disturb_nominal, traj) for traj in trajectories ]
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

    # Plot the Wind Direction Plots
    query_angles = collect(-π:0.01:π)

    polar_plot = plot(legend=false, framestyle = :box, dpi = 300, proj = :polar)

    # Plot the likelihood of the wind direction for the nominal
    wind_components_nominal = disturb_nominal.wind_dist.wind_dir_dists.components
    for (c_ind, component) in enumerate(wind_components_nominal)
        comp_likelihoods = wind_dir_likelihood_circ(
            component, query_angles)
        plot!(polar_plot, query_angles, comp_likelihoods, 
                lw=2, ls=:dot, color=c_ind)
    end

    # Plot the total likelihood
    total_likelihood_nominal = wind_dir_likelihood_circ(
        disturb_nominal.wind_dist, query_angles)
        plot!(polar_plot, query_angles, total_likelihood_nominal, 
                lw=2, color="gray")

    # Plot the likelihood of the wind direction for the fuzzed
    wind_components = disturb.wind_dist.wind_dir_dists.components
    for (c_ind, component) in enumerate(wind_components)
        comp_likelihoods = wind_dir_likelihood_circ(
            component, query_angles)
        plot!(polar_plot, query_angles, comp_likelihoods, 
                lw=2, color=c_ind)
    end

    # Plot the total likelihood
    total_likelihood = wind_dir_likelihood_circ(
        disturb.wind_dist, query_angles)
    plot!(polar_plot, query_angles, total_likelihood, 
            lw=2, color="black")

    savefig(polar_plot, 
        joinpath(save_dir, "wind_direction_likelihood_" * dataset * ".png"))

    println("=== Disturbance Example Completed ===\n")
end


# Instead of reloading the scene data, we now pass the BurnScene object.
# Note that the burn scene already includes the simulation number and 
# moisture level.
function run_burn_scene_with_background(burn_scene::BurnScene, dataset::String) #, simulation::Int, moisture_level::String)
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

# Instead of reloading the scene data, we now pass the BurnScene object.
# Note that the burn scene already includes the simulation number and 
# moisture level.
function run_changing_wind_scene(burn_scene::BurnScene, dataset::String) #, moisture_level::String)
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
    vert_offset = 0
    horizontal_offset = 150
    txt_size = 8

    num_timesteps = length(burn_scene.t) + 1

    # --- Plot Sensor Readings ---
    for t in 1:num_timesteps
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
            annotate!(pt[1] + horizontal_offset, pt[2] + vert_offset,
                    text(round(sensor_readings_per_time_log[t][i], digits=1), txt_size),
                    halign=:center, valign=:center, color="black")
        end
        
        xlims!(burn_scene.reference_bounds_x...)
        ylims!(burn_scene.reference_bounds_y...)
        curr_wind_dir = Int(wind_dirs[t])
        savefig(p, joinpath(sensor_save_dir, "example_sensor_readings_$(t)_$(curr_wind_dir).png"))
    end

    # --- Plot Perimeter Readings ---
    for t in 1:num_timesteps
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
            annotate!(pt[1] + horizontal_offset, pt[2] + vert_offset,
                    text(round(perimeter_readings_per_time_log[t][i], digits=1), txt_size),
                    halign=:center, valign=:center, color="black")
        end
        
        xlims!(burn_scene.reference_bounds_x...)
        ylims!(burn_scene.reference_bounds_y...)
        curr_wind_dir = Int(wind_dirs[t])
        savefig(p, joinpath(perimeter_save_dir, "example_perimeter_readings_$(t)_$(curr_wind_dir).png"))
    end

    # Plot the perimeter readings and sensor readings together as a 
    # side-by-side histograms per time step
    println("Plotting sensor vs perimeter histograms...")

    # lowest_loglikelihood = min(minimum(loglikelihoods_sampled), minimum(loglikelihoods_nominal))
    # highest_loglikelihood = max(maximum(loglikelihoods_sampled), maximum(loglikelihoods_nominal))
    # n_bins = 200
    # bin_edges = range(lowest_loglikelihood, stop=highest_loglikelihood, length=n_bins+1)

    # We want the x-axis to be the same for both histograms across all time steps
    # so we compute the bin edges once and use them for all histograms, but we
    # need to ignore the -Inf values in the sensor readings (due to log10)
    min_val_to_plot = -10.0
    max_val_to_plot = 0.0
    sensor_readings_all = vcat(sensor_readings_per_time_log...)
    # Print the number of -Inf values in the sensor readings
    println("Number of -Inf values in sensor readings: ", sum(sensor_readings_all .== -Inf))
    # Print the number of NaN values in the sensor readings
    println("Number of NaN values in sensor readings: ", sum(isnan.(sensor_readings_all)))
    sensor_readings_all = filter(x -> isfinite(x), sensor_readings_all)

    perimeter_readings_all = vcat(perimeter_readings_per_time_log...)
    # Print the number of -Inf values in the perimeter readings
    println("Number of -Inf values in perimeter readings: ", sum(perimeter_readings_all .== -Inf))
    # Print the number of NaN values in the perimeter readings
    println("Number of NaN values in perimeter readings: ", sum(isnan.(perimeter_readings_all)))
    perimeter_readings_all = filter(x -> isfinite(x), perimeter_readings_all)

    lowest_finite_reading = max(min_val_to_plot, 
        min(minimum(sensor_readings_all), minimum(perimeter_readings_all)))
    highest_finite_reading = max(max_val_to_plot,
        max(maximum(sensor_readings_all), maximum(perimeter_readings_all)))

    println("Lowest finite reading: ", lowest_finite_reading)
    println("Highest finite reading: ", highest_finite_reading)

    n_bins = 20
    bin_edges = range(lowest_finite_reading, stop=highest_finite_reading, length=n_bins+1)

    for t in 1:num_timesteps
        p_sensor = histogram(perimeter_readings_per_time_log[t], 
                            bins=bin_edges, label="Perimeter Readings",
                            xlabel="Log10 Perimeter Readings", ylabel="Frequency",
                            framestyle=:box, dpi=300, alpha=0.5, 
                            fillstyle = :/)
        histogram!(sensor_readings_per_time_log[t], 
                            bins=bin_edges, label="Sensor Readings",
                            xlabel="Log10 Sensor Readings", ylabel="Frequency",
                            framestyle=:box, dpi=300, alpha=0.5)

        # Check if failure (linear and log scales) and include it in title 
        # for clarity
        is_fail_linear = is_failure(sensor_readings_per_time[t],
                                    perimeter_readings_per_time[t],
                                    convert_to_linear=true)
        is_fail_log = is_failure(sensor_readings_per_time_log[t],
                                perimeter_readings_per_time_log[t],
                                convert_to_linear=false)
        fail_str_linear = is_fail_linear ? "FAIL (linear)" : "PASS (linear)"
        fail_str_log = is_fail_log ? "FAIL (log)" : "PASS (log)"

        title_str = "Time Step $t: $fail_str_linear, $fail_str_log"
        title!(title_str)
        
        savefig(p_sensor, joinpath(save_dir, "example_sensor_perimeter_histograms_$(t).png"))
    end
        

    
    # Plot plume maps over time
    println("Plotting plume maps...")
    bounds = vcat(burn_scene.reference_bounds_x..., burn_scene.reference_bounds_y...)
    for t in 1:num_timesteps
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



function run_most_likely_failure(burn_scene::BurnScene, 
        dataset::String, simulation::Int;
        num_trajectories::Int=100)
    println("\n=== Running Disturbance Example ===")

    # Reset Random seed for reproducibility
    Random.seed!(42)
    
    # Define the save directory using joinpath
    save_dir = joinpath("plots", "most_likely_failure", dataset)
    if !isdir(save_dir)
        mkpath(save_dir)
    end

    Q = burn_scene.Q

    # Choose disturbance distributions based on save_dir
    if dataset == "HenryCoe"
        disturb = HENRY_COE_FUZZED(Q)
        disturb_nominal = HENRY_COE_NOMINAL(Q)
    elseif dataset == "Malibu"
        disturb = MALIBU_FUZZED(Q)
        disturb_nominal = MALIBU_NOMINAL(Q)
    elseif dataset == "Shasta"
        disturb = SHASTA_FUZZED(Q)
        disturb_nominal = SHASTA_NOMINAL(Q)
    else
        println("Dataset not recognized. Using default disturbance.")
        disturb = DEFAULT_DISTURBANCE(Q)
        disturb_nominal = DEFAULT_DISTURBANCE(Q)
    end

    # Sample the distribution trajectories and compute likelihoods
    num_timesteps = length(burn_scene.t) + 1
    trajectories = [ sample_disturbances(disturb, num_timesteps) for _ in 1:num_trajectories ]
    loglikelihoods_sampled = [ disturbance_trajectory_log_likelihood(disturb, traj) for traj in trajectories ]
    loglikelihoods_nominal = [ disturbance_trajectory_log_likelihood(disturb_nominal, traj) for traj in trajectories ]

    # Get the sensor locations
    sensor_locations = burn_scene.snode_locations
    sensor_locations_3D = [ [pt[1], pt[2], 0.0] for pt in sensor_locations ]

    # Get the perimeter definitions
    perimeter_nodes = burn_scene.perimeter_sample_points
    perimeter_nodes_3D = [ [pt[1], pt[2], 0.0] for pt in perimeter_nodes ]

    # Define the failure arrays
    num_timesteps = length(burn_scene.t) + 1
    num_sensor_nodes = length(burn_scene.snode_locations)
    num_perimeter_nodes = length(perimeter_nodes)

    println("Number of timesteps: ", num_timesteps)
    println("Number of sensor nodes: ", num_sensor_nodes)
    println("Number of perimeter nodes: ", num_perimeter_nodes)

    sensor_readings = zeros(Float64, num_trajectories, num_sensor_nodes, num_timesteps)
    perimeter_readings = zeros(Float64, num_trajectories, num_perimeter_nodes, num_timesteps)

    is_fail_linear = zeros(Bool, num_trajectories, num_timesteps)
    is_fail_log = zeros(Bool, num_trajectories, num_timesteps)

    is_fail_linear_any = zeros(Bool, num_trajectories)
    is_fail_log_any = zeros(Bool, num_trajectories)
    is_fail_linear_all = zeros(Bool, num_trajectories)
    is_fail_log_all = zeros(Bool, num_trajectories)
    is_fail_linear_half = zeros(Bool, num_trajectories)
    is_fail_log_half = zeros(Bool, num_trajectories)

    num_time_fails_linear = zeros(Int, num_trajectories)
    num_time_fails_log = zeros(Int, num_trajectories)

    # Now run the simulations from the disturbances 
    for (traj_ind, traj) in enumerate(trajectories)
        h_plume = traj.height
        air_class = "D" # Could be variable in the future

        # Generate the smoke samples
        n_smoke_samples = 10
        smoke_samples_per_time = gen_smoke_samples_over_time(burn_scene, 
            n_smoke_samples=n_smoke_samples)

        # Generate the plumes
        plumes_per_time = gen_plumes_over_time(smoke_samples_per_time, 
            Q, h_plume, air_class)

        # Convert the wind directions and speeds to vectors
        wind_vecs = to_wind_vec(traj.wind_dirs, traj.wind_speeds)

        # Generate the sensor readings
        # Sensor readings are vectors of vectors, but we need to store them
        # in a matrix. So, we need to convert them to a matrix
        sensor_readings[traj_ind, :, :] = hcat(
            gen_sensor_readings_of_plume(
                plumes_per_time, sensor_locations_3D, wind_vecs)...)

        # Generate the perimeter readings
        perimeter_readings[traj_ind, :, :] = hcat(
            gen_sensor_readings_of_plume(
                plumes_per_time, perimeter_nodes_3D, wind_vecs)...)
    end

    # Convert the sensor readings to log10
    sensor_readings_log = log10.(sensor_readings)
    perimeter_readings_log = log10.(perimeter_readings)

    # Check if failure (linear and log scales)
    for traj_ind in 1:num_trajectories
        for t in 1:num_timesteps
            is_fail_linear[traj_ind, t] = is_failure(
                sensor_readings_log[traj_ind, :, t],
                perimeter_readings_log[traj_ind, :, t],
                convert_to_linear=true)
            is_fail_log[traj_ind, t] = is_failure(
                sensor_readings_log[traj_ind, :, t],
                perimeter_readings_log[traj_ind, :, t],
                convert_to_linear=false)
        end

        is_fail_linear_any[traj_ind] = any(is_fail_linear[traj_ind, :])
        is_fail_log_any[traj_ind] = any(is_fail_log[traj_ind, :])
        is_fail_linear_all[traj_ind] = all(is_fail_linear[traj_ind, :])
        is_fail_log_all[traj_ind] = all(is_fail_log[traj_ind, :])

        # Check if half of the timesteps fail
        is_fail_linear_half[traj_ind] = sum(is_fail_linear[traj_ind, :]) >= num_timesteps / 2
        is_fail_log_half[traj_ind] = sum(is_fail_log[traj_ind, :]) >= num_timesteps / 2

        num_time_fails_linear[traj_ind] = sum(is_fail_linear[traj_ind, :])
        num_time_fails_log[traj_ind] = sum(is_fail_log[traj_ind, :])
    end

    # Print the number of failures
    println("Number of trajectories: ", num_trajectories)
    println("Number of failures (linear, any): ", sum(is_fail_linear_any))
    println("Number of failures (log, any): ", sum(is_fail_log_any))
    println("Number of failures (linear, all): ", sum(is_fail_linear_all))
    println("Number of failures (log, all): ", sum(is_fail_log_all))
    println("Number of failures (linear, half): ", sum(is_fail_linear_half))
    println("Number of failures (log, half): ", sum(is_fail_log_half))

    println("Median number of timestep failures (linear): ", median(num_time_fails_linear))
    println("Median number of timestep failures (log): ", median(num_time_fails_log))
    println("Mean number of timestep failures (linear): ", mean(num_time_fails_linear))
    println("Mean number of timestep failures (log): ", mean(num_time_fails_log))

    # Save the above numbers to file
    save_file = joinpath(save_dir, "most_likely_failure_sim$(simulation).txt")
    open(save_file, "w") do io
        println(io, "Number of trajectories: ", num_trajectories)
        println(io, "Number of failures (linear, any): ", sum(is_fail_linear_any))
        println(io, "Number of failures (log, any): ", sum(is_fail_log_any))
        println(io, "Number of failures (linear, all): ", sum(is_fail_linear_all))
        println(io, "Number of failures (log, all): ", sum(is_fail_log_all))
        println(io, "Number of failures (linear, half): ", sum(is_fail_linear_half))
        println(io, "Number of failures (log, half): ", sum(is_fail_log_half))
        println(io, "Median number of timestep failures (linear): ", median(num_time_fails_linear))
        println(io, "Median number of timestep failures (log): ", median(num_time_fails_log))
        println(io, "Mean number of timestep failures (linear): ", mean(num_time_fails_linear))
        println(io, "Mean number of timestep failures (log): ", mean(num_time_fails_log))
    end

    # Filter the log-likelihoods based on the failure criteria
    # Sampled
    logL_fail_linear_any_sampled = loglikelihoods_sampled[is_fail_linear_any]
    logL_fail_log_any_sampled = loglikelihoods_sampled[is_fail_log_any]
    logL_fail_linear_all_sampled = loglikelihoods_sampled[is_fail_linear_all]
    logL_fail_log_all_sampled = loglikelihoods_sampled[is_fail_log_all]
    logL_fail_linear_half_sampled = loglikelihoods_sampled[is_fail_linear_half]
    logL_fail_log_half_sampled = loglikelihoods_sampled[is_fail_log_half]

    # Nominal
    logL_fail_linear_any_nominal = loglikelihoods_nominal[is_fail_linear_any]
    logL_fail_log_any_nominal = loglikelihoods_nominal[is_fail_log_any]
    logL_fail_linear_all_nominal = loglikelihoods_nominal[is_fail_linear_all]
    logL_fail_log_all_nominal = loglikelihoods_nominal[is_fail_log_all]
    logL_fail_linear_half_nominal = loglikelihoods_nominal[is_fail_linear_half]
    logL_fail_log_half_nominal = loglikelihoods_nominal[is_fail_log_half]

    # Get the top k failure likelihoods
    k_top = 10
    println("\n--- Top $k_top Disturbance Trajectories (Linear, Any) ---")
    if k_top < length(logL_fail_linear_any_nominal)
        println("\nLinear, Any:")
        top_indices_fail_linear_any_nominal = partialsortperm(
            logL_fail_linear_any_nominal, 1:k_top, rev=true)

        top_logL_fail_linear_any_nominal = logL_fail_linear_any_nominal[top_indices_fail_linear_any_nominal]
        top_logL_fail_linear_any_sampled = logL_fail_linear_any_sampled[top_indices_fail_linear_any_nominal]

        println("top indices: ", top_indices_fail_linear_any_nominal)
        println("log-likelihoods nominal: ", top_logL_fail_linear_any_nominal)
        println("log-likelihoods sampled: ", top_logL_fail_linear_any_sampled)
    else
        top_indices_fail_linear_any_nominal = ones(Int, k_top)
        top_logL_fail_linear_any_nominal = fill(NaN, k_top)
        top_logL_fail_linear_any_sampled = fill(NaN, k_top)
    end

    if k_top < length(logL_fail_log_any_nominal)
        println("\nLog, Any:")
        top_indices_fail_log_any_nominal = partialsortperm(
            logL_fail_log_any_nominal, 1:k_top, rev=true)

        top_logL_fail_log_any_nominal = logL_fail_log_any_nominal[top_indices_fail_log_any_nominal]
        top_logL_fail_log_any_sampled = logL_fail_log_any_sampled[top_indices_fail_log_any_nominal]

        println("top indices: ", top_indices_fail_log_any_nominal)
        println("log-likelihoods nominal: ", top_logL_fail_log_any_nominal)
        println("log-likelihoods sampled: ", top_logL_fail_log_any_sampled)
    else
        top_indices_fail_log_any_nominal = ones(Int, k_top)
        top_logL_fail_log_any_nominal = fill(NaN, k_top)
        top_logL_fail_log_any_sampled = fill(NaN, k_top)
    end

    if k_top < length(logL_fail_linear_all_nominal)
        println("\nLinear, All:")
        top_indices_fail_linear_all_nominal = partialsortperm(
            logL_fail_linear_all_nominal, 1:k_top, rev=true)

        top_logL_fail_linear_all_nominal = logL_fail_linear_all_nominal[top_indices_fail_linear_all_nominal]
        top_logL_fail_linear_all_sampled = logL_fail_linear_all_sampled[top_indices_fail_linear_all_nominal]

        println("top indices: ", top_indices_fail_linear_all_nominal)
        println("log-likelihoods nominal: ", top_logL_fail_linear_all_nominal)
        println("log-likelihoods sampled: ", top_logL_fail_linear_all_sampled)
    else
        top_indices_fail_linear_all_nominal = ones(Int, k_top)
        top_logL_fail_linear_all_nominal = fill(NaN, k_top)
        top_logL_fail_linear_all_sampled = fill(NaN, k_top)
    end

    if k_top < length(logL_fail_log_all_nominal)
        println("\nLog, All:")
        top_indices_fail_log_all_nominal = partialsortperm(
            logL_fail_log_all_nominal, 1:k_top, rev=true)

        top_logL_fail_log_all_nominal = logL_fail_log_all_nominal[top_indices_fail_log_all_nominal]
        top_logL_fail_log_all_sampled = logL_fail_log_all_sampled[top_indices_fail_log_all_nominal]

        println("top indices: ", top_indices_fail_log_all_nominal)
        println("log-likelihoods nominal: ", top_logL_fail_log_all_nominal)
        println("log-likelihoods sampled: ", top_logL_fail_log_all_sampled)
    else
        top_indices_fail_log_all_nominal = ones(Int, k_top)
        top_logL_fail_log_all_nominal = fill(NaN, k_top)
        top_logL_fail_log_all_sampled = fill(NaN, k_top)
    end

    if k_top < length(logL_fail_linear_half_nominal)
        println("\nLinear, Half:")
        top_indices_fail_linear_half_nominal = partialsortperm(
            logL_fail_linear_half_nominal, 1:k_top, rev=true)

        top_logL_fail_linear_half_nominal = logL_fail_linear_half_nominal[top_indices_fail_linear_half_nominal]
        top_logL_fail_linear_half_sampled = logL_fail_linear_half_sampled[top_indices_fail_linear_half_nominal]

        println("top indices: ", top_indices_fail_linear_half_nominal)
        println("log-likelihoods nominal: ", top_logL_fail_linear_half_nominal)
        println("log-likelihoods sampled: ", top_logL_fail_linear_half_sampled)
    else
        top_indices_fail_linear_half_nominal = ones(Int, k_top)
        top_logL_fail_linear_half_nominal = fill(NaN, k_top)
        top_logL_fail_linear_half_sampled = fill(NaN, k_top)
    end

    if k_top < length(logL_fail_log_half_nominal)
        println("\nLog, Half:")
        top_indices_fail_log_half_nominal = partialsortperm(
            logL_fail_log_half_nominal, 1:k_top, rev=true)

        top_logL_fail_log_half_nominal = logL_fail_log_half_nominal[top_indices_fail_log_half_nominal]
        top_logL_fail_log_half_sampled = logL_fail_log_half_sampled[top_indices_fail_log_half_nominal]

        # println(top_indices_fail_log_half)
        println("top indices: ", top_indices_fail_log_half_nominal)
        println("log-likelihoods nominal: ", top_logL_fail_log_half_nominal)
        println("log-likelihoods sampled: ", top_logL_fail_log_half_sampled)
    else
        top_indices_fail_log_half_nominal = ones(Int, k_top)
        top_logL_fail_log_half_nominal = fill(NaN, k_top)
        top_logL_fail_log_half_sampled = fill(NaN, k_top)
    end

    # Save to CSV the failure checks
    CSV.write(joinpath(save_dir, "failure_checks_sim$simulation.csv"), 
        DataFrame(
            traj_ind=1:num_trajectories,
            is_fail_linear_any=is_fail_linear_any,
            is_fail_log_any=is_fail_log_any,
            is_fail_linear_all=is_fail_linear_all,
            is_fail_log_all=is_fail_log_all,
            is_fail_linear_half=is_fail_linear_half,
            is_fail_log_half=is_fail_log_half,
            num_time_fails_linear=num_time_fails_linear,
            num_time_fails_log=num_time_fails_log
        )
    )

    # Save to CSV the likelihoods for the top k failures if they exist
    CSV.write(joinpath(save_dir, "top_k_failures_sim$simulation.csv"), 
        DataFrame(
            log_likelihood_fail_linear_any_nominal=top_logL_fail_linear_any_nominal,
            log_likelihood_fail_linear_any_sampled=top_logL_fail_linear_any_sampled,
            log_likelihood_fail_log_any_nominal=top_logL_fail_log_any_nominal,
            log_likelihood_fail_log_any_sampled=top_logL_fail_log_any_sampled,
            log_likelihood_fail_linear_all_nominal=top_logL_fail_linear_all_nominal,
            log_likelihood_fail_linear_all_sampled=top_logL_fail_linear_all_sampled,
            log_likelihood_fail_log_all_nominal=top_logL_fail_log_all_nominal,
            log_likelihood_fail_log_all_sampled=top_logL_fail_log_all_sampled,
            log_likelihood_fail_linear_half_nominal=top_logL_fail_linear_half_nominal,
            log_likelihood_fail_linear_half_sampled=top_logL_fail_linear_half_sampled,
            log_likelihood_fail_log_half_nominal=top_logL_fail_log_half_nominal,
            log_likelihood_fail_log_half_sampled=top_logL_fail_log_half_sampled
        )
    )


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

    # run_disturbance_example(burn_scene, dataset)
    # run_burn_scene_with_background(burn_scene, dataset) #, simulation, moisture_level)
    # run_changing_wind_scene(burn_scene, dataset) #, moisture_level)

    run_most_likely_failure(burn_scene, dataset, simulation)
end

main()

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
include("sample_perimeter.jl")

# Set a global RNG seed for reproducibility
Random.seed!(42)

###############################
# Helper Function to Compute Q
###############################
function compute_Q(dataset::String, moisture_level::String)
    Acreage_Consumption = joinpath("data", "FuelConsumptionProfiles", "$(dataset)_Scenarios.csv")
    df = CSV.read(Acreage_Consumption, DataFrame)
    consumed_total = parse(Float64, df[11, Symbol(moisture_level)])
    duration = 3600 * 24  # seconds in one day burn
    Q = consumed_total * 1e6 / duration
    println("Computed Q value for $(dataset): ", Q, " g/acre/sec")
    return Q
end

###############################
# Helper Functions for Scene and Simulation
###############################
function run_disturbance_example()
    println("\n=== Running Disturbance Example ===")
    
    # Define the save directory using joinpath
    save_dir = joinpath("plots", "disturbance_examples", "Malibu_fuzzed")
    if !isdir(save_dir)
        mkpath(save_dir)
    end

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
        HeightDistribution(400.0, 30.0),
        WindDistribution(
            [-120.0, 90.0, 150.0],
            [9.0, 5.0, 5.0],
            [0.1, 0.7, 0.2],
            10.0, 5.0;
            in_deg=true
        )
    )
    malibu_fuzzed_disturb = FullDisturbances(
        HeightDistribution(800.0, 5.0),
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

function run_burn_scene_with_background()
    println("\n=== Running Burn Scene with Background ===")
    # Model Input Parameters
    dataset = "Shasta"      # Change as needed
    simulation = 1          # 1-3
    moisture_level = "moist"  # "moderate", "moist", "wet"

    # Set up save directory using joinpath
    save_dir = joinpath("plots", dataset)
    if !isdir(save_dir)
        mkdir(save_dir)
    end

    # Load the background image
    img_path = joinpath("data", "BurnData", "DEMs_and_Buffered_Burns", "Background_$(dataset).png")
    background_img = load(img_path)
    background_img = RGB.(background_img)
    println("Background image dimensions: ", size(background_img))

    # Load burn scene data (without perimeter sampling)
    reference_bounds_filename = joinpath("data", "BurnData", dataset, "ReferencePoints", "referencelocations.txt")
    burnareas_filename = joinpath("data", "BurnData", dataset, "BurnAreas", "$(dataset)BurnAreas.txt")
    snode_locations_filename = joinpath("data", "BurnData", dataset, "SNodeLocations", "snodelocations_sim$(simulation).txt")
    burn_scene_obj = load_burn_scene_from_files(burnareas_filename, snode_locations_filename, reference_bounds_filename)

    println("Burn Scene Specs:")
    println(" - Burn polygons: ", length(burn_scene_obj.burn_polys_t))
    println(" - Time steps: ", length(burn_scene_obj.t))
    println(" - Sensor nodes: ", length(burn_scene_obj.snode_locations))
    println(" - X bounds: ", burn_scene_obj.reference_bounds_x)
    println(" - Y bounds: ", burn_scene_obj.reference_bounds_y)

    # Compute Q value using the helper function
    Q = compute_Q(dataset, moisture_level)

    # Create coordinate arrays for the background image
    xlims = burn_scene_obj.reference_bounds_x
    ylims = burn_scene_obj.reference_bounds_y
    nx = size(background_img, 2)
    ny = size(background_img, 1)
    background_x = range(xlims[1], stop=xlims[2], length=nx)
    background_y = range(ylims[1], stop=ylims[2], length=ny)

    # Plot the burn scene over time
    for t_ind in 1:(length(burn_scene_obj.t) + 1)
        p = plot_scene(burn_scene_obj, t_ind,
                       n_smoke_samples=10,
                       background_image=background_img,
                       background_x=background_x,
                       background_y=background_y)
        savefig(p, joinpath(save_dir, "burn_scene_test_$(dataset)_topobkgd_$(t_ind).png"))
    end

    println("=== Burn Scene with Background Completed ===\n")
end

function run_changing_wind_scene()
    println("\n=== Running Changing Wind Scene ===")
    dataset = "Malibu"
    moisture_level = "moist"  # Adjust if needed

    # Set up save directory using joinpath
    save_dir = joinpath("plots", dataset, "changing_wind")
    if !isdir(save_dir)
        mkdir(save_dir)
    end

    # Load burn scene data for Malibu with perimeter sampling.
    # Provide the total burn area filename so that perimeter_sample_points gets set.
    burnareas_filename = joinpath("data", "BurnData", dataset, "BurnAreas", "$(dataset)BurnAreas.txt")
    snode_locations_filename = joinpath("data", "BurnData", dataset, "SNodeLocations", "snodelocations_sim1.txt")
    ref_bounds_filename = joinpath("data", "BurnData", dataset, "ReferencePoints", "referencelocations.txt")
    total_burn_area_filename = joinpath("data", "BurnData", dataset, "BurnAreas","TotalBurnArea.txt")  # adjust as needed
    burn_scene_obj = load_burn_scene_from_files(burnareas_filename, snode_locations_filename, ref_bounds_filename, total_burn_area_filename)

    # Generate smoke samples and plumes
    n_smoke_samples = 10
    smoke_samples_per_time = gen_smoke_samples_over_time(burn_scene_obj, n_smoke_samples=n_smoke_samples)
    
    # Compute Q value for Malibu using the helper function
    Q = compute_Q(dataset, moisture_level)
    
    h_plume = 10.0
    air_class = "D"
    plumes_per_time = gen_plumes_over_time(smoke_samples_per_time, Q, h_plume, air_class)

    # Use the perimeter pseudo-nodes (sampled from the total burn area file) as sensor locations
    sensor_locations = burn_scene_obj.perimeter_sample_points
    # If sensor locations need to be extended to 3D, you could map:
    sensor_locations = [ [pt[1], pt[2], 0.0] for pt in sensor_locations ]

    # Compute sensor readings and then take log10
    sensor_readings_per_time = gen_sensor_readings_of_plume(plumes_per_time, sensor_locations, wind_vecs=to_wind_vec.( [30.0,35.0,40.0,35.0,-30.0,-25.0,-20.0,35.0,30.0] ))
    sensor_readings_per_time_log = [log10.(sensor_readings) for sensor_readings in sensor_readings_per_time]

    # For the wind vectors, ensure you define them (here we define them inline)
    wind_speed = 10.0
    wind_dirs = [30.0, 35.0, 40.0, 35.0, -30.0, -25.0, -20.0, 35.0, 30.0]
    to_wind_vec(wdir) = wind_speed .* [cosd(wdir), sind(wdir), 0.0]
    wind_vecs = to_wind_vec.(wind_dirs)
    wind_len_mult = 3.0
    vmin = -10.0
    vmax = 0.0

    # Plot sensor readings for each time step
    for t in 1:(length(burn_scene_obj.t) + 1)
        p = plot(framestyle=:box)
        curr_wind = wind_vecs[t]

        # Plot burn polygons with differentiated colors for past/current/future
        for ind in eachindex(burn_scene_obj.burn_polys_t)
            if ind < t
                plot!(burn_scene_obj.burn_polys_t[ind], color=:black, alpha=0.5, label="", linealpha=0.0)
            elseif ind == t
                plot!(burn_scene_obj.burn_polys_t[ind], color=:red, label="")
            else
                plot!(burn_scene_obj.burn_polys_t[ind], color=:lightgray, alpha=0.5, label="", linealpha=0.0)
            end
        end

        # Plot smoke sample locations and wind direction arrows
        if length(smoke_samples_per_time[t]) > 0
            for source_location in eachcol(smoke_samples_per_time[t])
                source_location = source_location[:]  # convert to 1D array
                scatter!([source_location[1]], [source_location[2]], color="gray", label="")
                quiver!([source_location[1]], [source_location[2]],
                        quiver=([wind_len_mult * curr_wind[1]], [wind_len_mult * curr_wind[2]]),
                        color="black", lw=2, label="")
            end
        end

        # Plot sensor readings with a gradient
        scatter!([pt[1] for pt in sensor_locations],
                 [pt[2] for pt in sensor_locations],
                 zcolor=[sensor_reading for sensor_reading in sensor_readings_per_time_log[t]],
                 label="",
                 c=cgrad(:bilbao, rev=true),
                 clim=(vmin, vmax))

        # Annotate sensor readings with values
        vert_offset = 30
        txt_size = 8
        for (i, pt) in enumerate(sensor_locations)
            annotate!(pt[1], pt[2] + vert_offset,
                      text(round(sensor_readings_per_time_log[t][i], digits=1), txt_size),
                      halign=:center, valign=:center, color="black")
        end

        xlims!(burn_scene_obj.reference_bounds_x...)
        ylims!(burn_scene_obj.reference_bounds_y...)
        curr_wind_dir = Int(wind_dirs[t])
        savefig(p, joinpath(save_dir, "example_readings_$(t)_$(curr_wind_dir).png"))
    end

    # Plot the plume maps over time
    println("Plotting plume maps...")
    bounds = vcat(burn_scene_obj.reference_bounds_x..., burn_scene_obj.reference_bounds_y...)
    for t in 1:(length(burn_scene_obj.t) + 1)
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
        xlims!(burn_scene_obj.reference_bounds_x...)
        ylims!(burn_scene_obj.reference_bounds_y...)
        curr_wind_dir = Int(wind_dirs[t])
        savefig(p, joinpath(save_dir, "example_plume_multiple_$(t)_$(curr_wind_dir).png"))
    end

    println("=== Changing Wind Scene Completed ===\n")
end

function main()
    run_disturbance_example()
    run_burn_scene_with_background()
    run_changing_wind_scene()
end

main()

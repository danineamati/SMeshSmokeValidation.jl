# Run the scene of the burn

using Plots
using Random

using SMeshSmokeValidation

# Set the RNG seed for reproducibility
Random.seed!(42)

dataset = "HenryCoe"

# Save directory for the plots
save_dir = "plots/" * dataset * "/changing_wind"
if !isdir(save_dir)
    mkdir(save_dir)
end

# Parse the burn areas
burnareas_filename = "data\\BurnData\\" * dataset * "\\BurnAreas\\" * dataset * "BurnAreas.txt"
snode_locations_filename = "data\\BurnData\\" * dataset * "\\SNodeLocations\\snodelocations.txt"
ref_bounds_filename = "data\\BurnData\\" * dataset * "\\ReferencePoints\\referencelocations.txt"

burn_scene_obj = load_burn_scene_from_files(
    burnareas_filename, snode_locations_filename, ref_bounds_filename)

# Generate the smoke samples
n_smoke_samples = 10
smoke_samples_per_time = gen_smoke_samples_over_time(burn_scene_obj, 
        n_smoke_samples=n_smoke_samples)

# Generate the plumes
Q = 1.0
h_plume = 10.0
air_class = "D"
plumes_per_time = gen_plumes_over_time(smoke_samples_per_time, 
                        Q, h_plume, air_class)

# Get the sensor readings
wind_speed = 10.0
wind_dirs = [30.0, 35.0, 40.0, 35.0, -30.0, -25.0, -20.0]
to_wind_vec(wdir) = wind_speed .* [cosd(wdir), sind(wdir), 0.0]
wind_vecs = to_wind_vec.(wind_dirs)

wind_len_mult = 3.0

vmin = -10.0
vmax = 0.0

# Sensor locations need to be 3D
sensor_locations = burn_scene_obj.snode_locations
for s_ind in eachindex(sensor_locations)
    sensor_locations[s_ind] = [sensor_locations[s_ind]; 0.0]
end


sensor_readings_per_time = gen_sensor_readings_of_plume(plumes_per_time, 
    burn_scene_obj.snode_locations, wind_vecs)

sensor_readings_per_time_log = [log10.(sensor_readings_at_time)
    for sensor_readings_at_time in sensor_readings_per_time]

for t in 1:(length(burn_scene_obj.t) + 1)

    p = plot(framestyle = :box)
    curr_wind = wind_vecs[t]

    # Plot the polygons in order
    for ind in eachindex(burn_scene_obj.burn_polys_t)
        if ind < t
            plot!(burn_scene_obj.burn_polys_t[ind], color=:black, alpha=0.5, label="",
                    linealpha=0.0)
        elseif ind == t
            plot!(burn_scene_obj.burn_polys_t[ind], color=:red, label="")
        else
            plot!(burn_scene_obj.burn_polys_t[ind], color=:lightgray, alpha=0.5, label="",
                    linealpha=0.0)
        end
    end

    if length(smoke_samples_per_time[t]) > 0
        for source_location in eachcol(smoke_samples_per_time[t])
            # Cast to a 1D array
            source_location = source_location[:]

            # Add a scatter for the source location
            scatter!([source_location[1]], [source_location[2]], 
                    color="gray", label="")

            # Add a quiver for the wind direction
            quiver!([source_location[1]], [source_location[2]], 
                    quiver=([wind_len_mult * curr_wind[1]], [wind_len_mult * curr_wind[2]]),
                    color="black", lw=2, label="")        
        end
    end

    # Add the sensor readings
    scatter!([sensor[1] for sensor in sensor_locations], 
        [sensor[2] for sensor in sensor_locations], 
        zcolor=[sensor_reading for sensor_reading in sensor_readings_per_time_log[t]], 
        label="",
        c=cgrad(:bilbao, rev=true),
        clim=(vmin, vmax))

    # Add text for the sensor readings
    vert_offset = 30
    txt_size = 8
    for (i, sensor_loc) in enumerate(sensor_locations)
        annotate!(sensor_loc[1], sensor_loc[2] + vert_offset, 
            text(round(sensor_readings_per_time_log[t][i], digits=1), txt_size), 
            halign=:center, valign=:center, color="black")
    end

    # Enforce the plot area bounds
    xlims!(burn_scene_obj.reference_bounds_x...)
    ylims!(burn_scene_obj.reference_bounds_y...)

    curr_wind_dir = Int(wind_dirs[t])
    savefig(p, joinpath(save_dir, "example_readings_$(t)_$(curr_wind_dir).png"))
end

# Repeat, but now plot the plumes
println("Plotting plumes...")

bounds = vcat(burn_scene_obj.reference_bounds_x..., burn_scene_obj.reference_bounds_y...)

for t in 1:(length(burn_scene_obj.t) + 1)

    curr_wind = wind_vecs[t]

    p = plot_multiple_plumes_bounds(plumes_per_time[t], bounds, curr_wind, 
        x_num=101, y_num=103, vmax=vmax)

    if length(smoke_samples_per_time[t]) > 0
        for source_location in eachcol(smoke_samples_per_time[t])
            # Cast to a 1D array
            source_location = source_location[:]

            # Add a scatter for the source location
            scatter!([source_location[1]], [source_location[2]], 
                    color="gray", label="")

            # Add a quiver for the wind direction
            quiver!([source_location[1]], [source_location[2]], 
                    quiver=([wind_len_mult * curr_wind[1]], [wind_len_mult * curr_wind[2]]),
                    color="black", lw=2, label="")        
        end
    end

    # Enforce the plot area bounds
    xlims!(burn_scene_obj.reference_bounds_x...)
    ylims!(burn_scene_obj.reference_bounds_y...)

    curr_wind_dir = Int(wind_dirs[t])
    savefig(p, joinpath(save_dir, "example_plume_multiple_$(t)_$(curr_wind_dir).png"))
    
end

println("Done!")


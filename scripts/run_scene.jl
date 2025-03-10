# Run the scene of the burn

using SMeshSmokeValidation

###########################################
# Functions to move to library later
###########################################


function gen_smoke_samples_over_time(burn_scene_obj::BurnScene;
    n_smoke_samples::Int=10)
    # Extract the polygons and time values from the burn_scene object
    burn_polys_t = burn_scene_obj.burn_polys_t

    smoke_samples_per_time = []
    for t_ind in 1:(1 + length(burn_polys_t))
        smoke_samples = sample_smoke_from_scene(burn_scene_obj, n_smoke_samples,
                                                t_ind)
        smoke_samples_mat = hcat(smoke_samples...)
        push!(smoke_samples_per_time, smoke_samples_mat)
    end
    return smoke_samples_per_time
end


function gen_plumes_at_time(smoke_samples_at_time, Q, h_plume, air_class)
    plumes = Vector{PlumeFromPointSource}()

    for smoke_sample in eachcol(smoke_samples_at_time)
        # Smoke sample should be a 3D point
        # Cast to a 1D array
        smoke_sample = smoke_sample[:]

        # If 2D, convert to 3D
        if length(smoke_sample) == 2
            smoke_sample = [smoke_sample; 0.0]
        end

        plume = PlumeFromPointSource(Q, h_plume, air_class, smoke_sample)
        push!(plumes, plume)
    end
    return plumes
end


function gen_plumes_over_time(smoke_samples_per_time, Q, h_plume, air_class)
    plumes_per_time = Vector{Vector{PlumeFromPointSource}}()
    for smoke_samples_at_time in smoke_samples_per_time
        # If there are no smoke samples, then there are no plumes
        if isempty(smoke_samples_at_time)
            push!(plumes_per_time, [PlumeFromPointSource()])
            continue
        end

        plumes = gen_plumes_at_time(smoke_samples_at_time, Q, h_plume, air_class)
        push!(plumes_per_time, plumes)
    end
    return plumes_per_time
end


function gen_sensor_readings_of_plume(plumes_per_time, sensor_locations, wind_vec)
    sensor_readings_per_time = []
    for plumes_at_curr_time in plumes_per_time
        # If there are no plumes, then there are no sensor readings
        if isempty(plumes_at_curr_time)
            push!(sensor_readings_per_time, [])
            continue
        end

        sensor_readings = query_multiple_plumes_points(
            plumes_at_curr_time, sensor_locations, wind_vec)

        push!(sensor_readings_per_time, sensor_readings)

    end

    return sensor_readings_per_time
    
end


###########################################
# Script
###########################################

using Plots
using Random

# Set the RNG seed for reproducibility
Random.seed!(42)

dataset = "HenryCoe"

# Save directory for the plots
save_dir = "plots/" * dataset * "/test_snode_samples"
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
Q = 0.1
h_plume = 10.0
air_class = "D"
plumes_per_time = gen_plumes_over_time(smoke_samples_per_time, 
                        Q, h_plume, air_class)

# Get the sensor readings
wind_speed = 10.0
wind_dir = 30.0
wind_vec = wind_speed * [cosd(wind_dir), sind(wind_dir), 0.0]
wind_len_mult = 3.0

vmin = -10.0
vmax = 0.0

# Sensor locations need to be 3D
sensor_locations = burn_scene_obj.snode_locations
for s_ind in eachindex(sensor_locations)
    sensor_locations[s_ind] = [sensor_locations[s_ind]; 0.0]
end


sensor_readings_per_time = gen_sensor_readings_of_plume(plumes_per_time, 
    burn_scene_obj.snode_locations, wind_vec)

sensor_readings_per_time_log = [log10.(sensor_readings_at_time)
    for sensor_readings_at_time in sensor_readings_per_time]

for t in 1:(length(burn_scene_obj.t) + 1)

    p = plot()

    if length(smoke_samples_per_time[t]) > 0
        for source_location in eachcol(smoke_samples_per_time[t])
            # Cast to a 1D array
            source_location = source_location[:]

            # Add a scatter for the source location
            scatter!([source_location[1]], [source_location[2]], 
                    color="gray", label="")

            # Add a quiver for the wind direction
            quiver!([source_location[1]], [source_location[2]], 
                    quiver=([wind_len_mult * wind_vec[1]], [wind_len_mult * wind_vec[2]]),
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

    savefig(p, joinpath(save_dir, "example_plume_multiple_$(wind_dir)_$(t).png"))
end

println("Done!")


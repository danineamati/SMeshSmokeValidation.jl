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
    plumes = []
    for smoke_samples in smoke_samples_at_time
        plume = PlumeFromPointSource(Q, h_plume, air_class, smoke_samples)
        push!(plumes, plume)
    end
    return plumes
end


function gen_plumes_over_time(smoke_samples_per_time, Q, h_plume, air_class)
    plumes_per_time = []
    for smoke_samples_at_time in smoke_samples_per_time
        plumes = gen_plumes_at_time(smoke_samples_at_time, Q, h_plume, air_class)
        push!(plumes_per_time, plumes)
    end
    return plumes_per_time
end


function gen_sensor_readings_of_plume(plumes_per_time, sensor_locations, wind_vec)
    sensor_readings_per_time = []
    for plumes_at_curr_time in plumes_per_time

        sensor_readings = query_multiple_plumes_points(
            plumes_at_curr_time, sensor_locations, wind_vec)

        push!(sensor_readings_per_time, sensor_readings)

    end

    return sensor_readings_per_time
    
end


###########################################
# Script
###########################################

using Random


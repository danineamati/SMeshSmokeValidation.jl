
using Plots

using SMeshSmokeValidation

# Make plotting directory
save_dir = "plots/example_plume"
if !isdir(save_dir)
    mkdir(save_dir)
end

# Consider an example sensor in some example wind, check that we can indeed
# get the plume model over some bounds.

# Make plume model
Q = 100.0 # g/s
h_plume = 50.0 # m
wind_speed = 10.0 # m/s
wind_dirs = [0.0, 150.0, -30.0] # degrees
air_class = "C"
source_locations = [
    [300.0, -200.0, -50.0],
    [-100.0, 0.0, -20.0],
    [-120, 210.0, 50.0]
]
vmax = 0.0

# Make bounds
bounds = [-700.0, 700.0, -500.0, 500.0]

println("\nStarting the single plume plotting routine...\n")

for wind_dir in wind_dirs
    wind_vec = wind_speed * [cosd(wind_dir), sind(wind_dir), 0.0]


    for source_location in source_locations

        plume_ex = PlumeFromPointSource(Q, h_plume, air_class, source_location)

        # Plot the plume
        p = plot_single_plume_bounds(plume_ex, bounds, wind_vec, 
                                     x_num=101, y_num=103, vmax=vmax)

        # Add a scatter for the source location
        scatter!([source_location[1]], [source_location[2]], 
                color="red", label="")

        # Add a quiver for the wind direction
        quiver!([source_location[1]], [source_location[2]], 
                quiver=([10 * wind_vec[1]], [10 * wind_vec[2]]),
                color="black", lw=2, label="")

        # Save the plot
        savefig(p, joinpath(save_dir, "example_plume_$(source_location)_$(wind_dir).png"))
    end
end



println("\nStarting the multiple plume plotting routine...\n")


plume_vector = [PlumeFromPointSource(Q, h_plume, air_class, source_location) 
                    for source_location in source_locations]


for wind_dir in wind_dirs
    wind_vec = wind_speed * [cosd(wind_dir), sind(wind_dir), 0.0]

    # Plot the plumes
    p = plot_multiple_plumes_bounds(plume_vector, bounds, wind_vec, 
                                     x_num=101, y_num=103, vmax=vmax)

    for source_location in source_locations
        # Add a scatter for the source location
        scatter!([source_location[1]], [source_location[2]], 
                color="red", label="")

        # Add a quiver for the wind direction
        quiver!([source_location[1]], [source_location[2]], 
                quiver=([10 * wind_vec[1]], [10 * wind_vec[2]]),
                color="black", lw=2, label="")        
    end

    # Save the plot
    savefig(p, joinpath(save_dir, "example_plume_multiple_$(wind_dir).png"))
end

println("\nPlotting plumes with SNodes in a circle...\n")
# Make a circle of sensors
n_sensors = 20
r_sensors = 200.0
sensor_locations = [
    [r_sensors * cosd(360 * i / n_sensors), r_sensors * sind(360 * i / n_sensors), 0.0]
    for i in 1:n_sensors
]

for wind_dir in wind_dirs
    wind_vec = wind_speed * [cosd(wind_dir), sind(wind_dir), 0.0]

    # Plot the plumes
    # p = plot_multiple_plumes_bounds(plume_vector, bounds, wind_vec, 
    #                                  x_num=101, y_num=103, vmax=vmax)
    p = plot()

    for source_location in source_locations
        # Add a scatter for the source location
        scatter!([source_location[1]], [source_location[2]], 
                color="red", label="")

        # Add a quiver for the wind direction
        quiver!([source_location[1]], [source_location[2]], 
                quiver=([10 * wind_vec[1]], [10 * wind_vec[2]]),
                color="black", lw=2, label="")        
    end

    println("\nQuerying the plume model at the sensor locations...")

    # Get the sensor readings
    vmin = -10.0
    sensor_readings_raw = query_multiple_plumes_points(
        plume_vector, sensor_locations, wind_vec)

    sensor_readings_raw[sensor_readings_raw .<= 10^vmin] .= 10^vmin
    sensor_readings = log10.(sensor_readings_raw)

    # Add the sensor readings
    scatter!([sensor[1] for sensor in sensor_locations], 
             [sensor[2] for sensor in sensor_locations], 
             zcolor=[sensor_reading for sensor_reading in sensor_readings], 
             label="",
             c=cgrad(:bilbao, rev=true),
             clim=(vmin, vmax))

    # Add text for the sensor readings
    vert_offset = 20
    txt_size = 8
    for i in 1:n_sensors
        annotate!(sensor_locations[i][1], sensor_locations[i][2] + vert_offset, 
              text(round(sensor_readings[i], digits=1), txt_size), 
              halign=:center, valign=:center, color="black")
    end

    # Save the plot
    savefig(p, joinpath(save_dir, "example_plume_multiple_sensors_$(wind_dir).png"))


    # Plot a histogram of the sensor readings
    # Bins to vmin to vmax in 0.5 increments
    bins = collect(vmin:0.5:vmax)
    p_box = histogram(sensor_readings, label="Sensor Readings", 
                    xlabel="Sensor Number", ylabel="Log10 Sensor Reading",
                    title="Sensor Readings for Wind Direction $(wind_dir)",
                    bins=bins)

    savefig(p_box, joinpath(save_dir, "example_plume_sensors_histogram_$(wind_dir).png"))

end






println("Done plotting plumes.")

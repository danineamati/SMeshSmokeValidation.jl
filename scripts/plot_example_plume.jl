
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

println("Starting the single plume plotting routine...")

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



println("Starting the multiple plume plotting routine...")


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

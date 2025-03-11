# Define, sample, and plot the wind distribution
# 
# Mixture of distributions (i.e., mixture of VonMises Distributions) for the 
# direction and a normal distribution for the speed

using Distributions
using Plots

save_dir = "plots/wind_distribution_testing"
if !isdir(save_dir)
    mkdir(save_dir)
end


mixture_components = [
    VonMises(0.0, 9.0),
    VonMises(deg2rad(-120), 5.0),
    VonMises(deg2rad(-180), 5.0)
]

# Print the variance of each component
for (i, component) in enumerate(mixture_components)
    println("Component ", i, " standard deviation (deg): ", rad2deg(std(component)))
end

wind_dir_dists = MixtureModel(
    mixture_components,
    [0.7, 0.15, 0.15]
)

# Visualize the likelihood of the wind direction on a polar plot
polar_plot = plot(legend=false, framestyle = :box, proj = :polar)

# Plot the likelihood of the wind direction
# The Distributions.jl implementation of the VonMises distribution does not
# wrap around the circle, so we need to plot the likelihood for double the 
# circle and then collapse it to the circle (i.e., sum the likelihoods for
# the same angle but with branches from -2π to 0 and 0 to 2π)
query_angles = [collect(-2π:0.01:0), collect(0:0.01:2π)]
for component in mixture_components
    branch1_pdf = pdf.(component, query_angles[1])
    branch2_pdf = pdf.(component, query_angles[2])
    pdf_sum = branch1_pdf + branch2_pdf
    plot!(query_angles[2], pdf_sum, lw=2)
end

branch1_pdf = pdf.(wind_dir_dists, query_angles[1])
branch2_pdf = pdf.(wind_dir_dists, query_angles[2])
pdf_sum = branch1_pdf + branch2_pdf
plot!(polar_plot, query_angles[2], pdf_sum, lw=2, color="black")

savefig(polar_plot, joinpath(save_dir, "wind_direction_likelihood.png"))



# Include the wind speed

wind_speed_dist = Normal(5.0, 1.0)


# Visualize through sampling
n_samples = 1000


# Sample the wind directions
wind_dirs = rand(wind_dir_dists, n_samples)

# Sample the wind speeds
wind_speeds = rand(wind_speed_dist, n_samples)

# Convert the wind directions and speeds to vectors
to_wind_vec(wdir, wspeed) = wspeed .* [cos(wdir), sin(wdir)]

wind_vecs = to_wind_vec.(wind_dirs, wind_speeds)

# Plot the wind vectors
p = plot(framestyle = :box, proj = :polar)

for (wd, ws) in zip(wind_dirs, wind_speeds)
    # quiver!([0.0], [0.0], quiver=([wind_vec[1]], [wind_vec[2]]), 
    #         color="black", lw=2, label="", alpha=25 / n_samples)
    quiver!([0.0], [0.0], quiver=([wd], [ws]), 
            color="black", lw=2, label="", alpha=25 / n_samples)
end

# plot!(p, ratio=1, legend=false, xlims=(-10, 10), ylims=(-10, 10), framestyle = :box, grid=false, size=(800, 800))

savefig(p, joinpath(save_dir, "wind_distribution.png"))


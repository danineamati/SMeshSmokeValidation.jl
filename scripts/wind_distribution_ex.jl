# Define, sample, and plot the wind distribution
# 
# Mixture of distributions (i.e., mixture of VonMises Distributions) for the 
# direction and a normal distribution for the speed

using Distributions
using Plots
using Random

include("../src/wind_distribution.jl")

Random.seed!(42)

# save_dir = "plots/wind_distribution_testing"
save_dir = "plots/wind_examples/Malibu_nominal"
if !isdir(save_dir)
    mkpath(save_dir)
end


default_wind_dist = WindDistribution(
    [0.0, -120, -180], 
    [9.0, 5.0, 5.0], 
    [0.7, 0.15, 0.15], 
    5.0, 1.0;
    in_deg=true
)

# Malibu Wind Nominal
malibu_nominal_wind_dist = WindDistribution(
    [-120.0, # Towards the ocean following a Santa Ana event 
       90.0, # From the ocean
      150.0], # Along the slopes
    [9.0, 5.0, 5.0], 
    [0.1, 0.7, 0.2], 
    10.0, 5.0;
    in_deg=true
)

# Malibu Wind Fuzzed
malibu_fuzzed_wind_dist = WindDistribution(
    [-120.0, # Towards the ocean following a Santa Ana event 
       90.0, # From the ocean
      150.0], # Along the slopes
    [50.0, 5.0, 5.0], 
    [0.7, 0.25, 0.05], 
    15.0, 1.0;
    in_deg=true
)



if save_dir == "plots/wind_distribution_testing"
    wind_dist = default_wind_dist
    wind_dist_nominal = default_wind_dist

elseif save_dir == "plots/wind_examples/Malibu_nominal"
    wind_dist = malibu_nominal_wind_dist
    wind_dist_nominal = malibu_nominal_wind_dist

elseif save_dir == "plots/wind_examples/Malibu_fuzzed"
    wind_dist = malibu_fuzzed_wind_dist
    wind_dist_nominal = malibu_nominal_wind_dist
end


println("Wind Distribution: ", wind_dist)
println("Type of wind_dist: ", typeof(wind_dist))
println("Type of wind_dist.wind_dir_dists: ", typeof(wind_dist.wind_dir_dists))
println("Type of wind_dist.wind_speed_dist: ", typeof(wind_dist.wind_speed_dist))

# Print the standard deviation of each component
println("Standard deviations of the components [deg]: ", dir_stds(wind_dist))

# Visualize the likelihood of the wind direction on a polar plot
plot_types = [:cartesian, :polar]
plot_yscales = [:linear, :log]
query_angles = collect(-π:0.01:π)

for plot_type in plot_types
    for yscale in plot_yscales
        polar_plot = plot(legend=false, framestyle = :box, dpi = 300,
                        proj = plot_type, yscale=yscale)

        # Plot the likelihood of the wind direction
        for component in wind_dist.wind_dir_dists.components
            comp_likelihoods = wind_dir_likelihood_circ(
                component, query_angles)
            plot!(polar_plot, query_angles, comp_likelihoods, lw=2)

            # println("length of comp_likelihoods: ", length(comp_likelihoods))
            # println("comp_likelihoods between ", minimum(comp_likelihoods), 
            #     " and ", maximum(comp_likelihoods))
        end

        total_likelihood = wind_dir_likelihood_circ(wind_dist, query_angles)
        plot!(polar_plot, query_angles, total_likelihood, lw=2, color="black")

        # println("total_likelihood between ", minimum(total_likelihood), 
        #     " and ", maximum(total_likelihood))

        savefig(polar_plot, joinpath(save_dir, "wind_direction_likelihood_" * string(plot_type) * "_yscale_" * string(yscale) * ".png"))
    end
end

# Repeat for the log likelihood of the wind direction

log_min = Inf

for plot_type in plot_types
    polar_plot_log = plot(legend=false, framestyle = :box, dpi = 300, 
                            proj = plot_type)

    # Offset the log likelihood to avoid negative values
    if plot_type == :cartesian
        log_offset = 0.0
    else
        log_offset = log_min
    end
    println("log_offset: ", log_offset)

    for component in wind_dist.wind_dir_dists.components
        comp_loglikelihood = wind_dir_loglikelihood_circ(
            component, query_angles) .- log_offset

        global log_min = min(log_min, minimum(comp_loglikelihood))

        plot!(polar_plot_log, query_angles, comp_loglikelihood, lw=2)
    end

    total_loglikelihood = wind_dir_loglikelihood_circ(wind_dist, query_angles)
    plot!(polar_plot_log, query_angles, total_loglikelihood .- log_offset, 
            lw=2, color="black")

    savefig(polar_plot_log, joinpath(save_dir, "wind_direction_loglikelihood_" * string(plot_type) * ".png"))
end


# Visualize through sampling
n_samples = 1000

# Sample the wind trajectory
wind_dir_traj, wind_speed_traj = sample_wind_dir_speed_trajectory(
    wind_dist, n_samples)

# Plot the wind vectors
p = plot(framestyle = :box, proj = :polar, dpi = 300)

for (wd, ws) in zip(wind_dir_traj, wind_speed_traj)
    # quiver!([0.0], [0.0], quiver=([wind_vec[1]], [wind_vec[2]]), 
    #         color="black", lw=2, label="", alpha=25 / n_samples)
    quiver!([0.0], [0.0], quiver=([wd], [ws]), 
            color="black", lw=2, label="", alpha=25 / n_samples)
end

# plot!(p, ratio=1, legend=false, xlims=(-10, 10), ylims=(-10, 10), framestyle = :box, grid=false, size=(800, 800))

savefig(p, joinpath(save_dir, "wind_distribution.png"))


# Get the log likelihood of the wind trajectory
loglikelihoods_from_sampled = [wind_loglikelihood(wind_dist, wd, ws) 
    for (wd, ws) in zip(wind_dir_traj, wind_speed_traj)]

loglikelihoods_from_nominal = [wind_loglikelihood(wind_dist_nominal, wd, ws) 
    for (wd, ws) in zip(wind_dir_traj, wind_speed_traj)]

# Plot the log likelihoods as a histogram
lowest_loglikelihood = min(minimum(loglikelihoods_from_sampled), minimum(loglikelihoods_from_nominal))
highest_loglikelihood = max(maximum(loglikelihoods_from_sampled), maximum(loglikelihoods_from_nominal))

n_bins = 200
bin_edges = range(lowest_loglikelihood, stop=highest_loglikelihood, length=n_bins+1)


sampled_label = occursin("fuzzed", lowercase(save_dir)) ? "Fuzzed" : "Sampled"

h = histogram(loglikelihoods_from_nominal, label="Nominal",
    alpha=0.5, bins=bin_edges, normalize=:pdf, 
    xlabel="Wind Vector Log Likelihood", ylabel="Estimated PDF", 
    framestyle = :box, dpi = 300)
histogram!(loglikelihoods_from_sampled, label=sampled_label, 
    alpha=0.5, bins=bin_edges, normalize=:pdf)

savefig(h, joinpath(save_dir, "wind_loglikelihood_histogram.png"))

println("Done!")

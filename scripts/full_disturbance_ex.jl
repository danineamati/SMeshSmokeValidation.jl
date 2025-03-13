
using Plots
using Random

include("../src/wind_distribution.jl")
include("../src/height_distribution.jl")
include("../src/disturbance.jl")


Random.seed!(42)


save_dir = "plots/disturbance_examples/Malibu_fuzzed"
if !isdir(save_dir)
    mkpath(save_dir)
end

# Defaults
default_disturb = FullDisturbances(
    HeightDistribution(100.0, 80.0),
    WindDistribution(
        [0.0, -120, -180], 
        [9.0, 5.0, 5.0], 
        [0.7, 0.15, 0.15], 
        5.0, 1.0;
        in_deg=true
    )
)


# Malibu Nominal
malibu_nominal_disturb = FullDisturbances(
    HeightDistribution(40.0, 30.0),
    WindDistribution(
        [-120.0, 90.0, 150.0], 
        [9.0, 5.0, 5.0], 
        [0.1, 0.7, 0.2], 
        10.0, 5.0;
        in_deg=true
    )
)


# Malibu Wind Fuzzed
malibu_fuzzed_disturb = FullDisturbances(
    HeightDistribution(80.0, 5.0),
    WindDistribution(
        [-120.0, 90.0, 150.0], 
        [9.0, 5.0, 5.0], 
        [0.7, 0.15, 0.15], 
        5.0, 1.0;
        in_deg=true
    )
)


if save_dir == "plots/disturbance_examples/full_disturbance_testing"
    disturb = default_disturb
    disturb_nominal = default_disturb
elseif save_dir == "plots/disturbance_examples/Malibu_nominal"
    disturb = malibu_nominal_disturb
    disturb_nominal = malibu_nominal_disturb
elseif save_dir == "plots/disturbance_examples/Malibu_fuzzed"
    disturb = malibu_fuzzed_disturb
    disturb_nominal = malibu_nominal_disturb
end


# Sample a bundle of trajectories of disturbances and
# calculate the likelihood and log-likelihood of each
# trajectory
num_timesteps = 10
num_trajectories = 1000


# Trajectories are iid
trajectories = [
    sample_disturbances(disturb, num_timesteps) for _ in 1:num_trajectories
]

# Likelihoods
likelihoods_sampled = [
    disturbance_trajectory_likelihood(disturb, traj) for traj in trajectories
]
likelihoods_nominal = [
    disturbance_trajectory_likelihood(disturb_nominal, traj) for traj in trajectories
]

# Log-likelihoods
loglikelihoods_sampled = [
    disturbance_trajectory_log_likelihood(disturb, traj) for traj in trajectories
]
loglikelihoods_nominal = [
    disturbance_trajectory_log_likelihood(disturb_nominal, traj) for traj in trajectories
]


# Plot the log-likelihoods as histograms
lowest_loglikelihood = min(minimum(loglikelihoods_sampled), minimum(loglikelihoods_nominal))
highest_loglikelihood = max(maximum(loglikelihoods_sampled), maximum(loglikelihoods_nominal))


n_bins = 200
bin_edges = range(lowest_loglikelihood, stop=highest_loglikelihood, length=n_bins+1)

sampled_label = occursin("fuzzed", lowercase(save_dir)) ? "Fuzzed" : "Sampled"

h = histogram(loglikelihoods_nominal, label="Nominal",
    alpha=0.5, bins=bin_edges, normalize=:pdf, 
    xlabel="Disturbance Trajectory Log Likelihood", ylabel="Estimated PDF", 
    framestyle = :box, dpi = 300)
histogram!(loglikelihoods_sampled, label=sampled_label, 
    alpha=0.5, bins=bin_edges, normalize=:pdf)

savefig(h, joinpath(save_dir, "disturbance_loglikelihood_histogram.png"))


# Print the "k" disturbance trajectory with the highest log-likelihoods
# under the nominal distribution
k_top = 10
top_indices = sortperm(loglikelihoods_nominal, rev=true)[1:k_top]

println("Top $k_top disturbance trajectories under the nominal distribution:")
for i in 1:k_top
    println("\nTrajectory $i has log-likelihood ", loglikelihoods_nominal[top_indices[i]])
    println("Trajectory $i: ", trajectories[top_indices[i]])
end


println("Done!")


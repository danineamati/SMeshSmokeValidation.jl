# Define, sample, and plot the height distribution


using Distributions
using Plots
using Random

include("../src/height_distribution.jl")

Random.seed!(42)

save_dir = "plots/height_examples/Malibu_fuzzed"
if !isdir(save_dir)
    mkpath(save_dir)
end


default_height_dist = HeightDistribution(1000.0, 80.0)

# Malibu Height Nominal (moist)
malibu_nominal_height_dist = HeightDistribution(400.0, 30.0)

# Malibu Height Fuzzed (dry)
malibu_fuzzed_height_dist = HeightDistribution(800.0, 5.0)


if save_dir == "plots/height_examples/height_distribution_testing"
    height_dist = default_height_dist
    height_dist_nominal = default_height_dist

elseif save_dir == "plots/height_examples/Malibu_nominal"
    height_dist = malibu_nominal_height_dist
    height_dist_nominal = malibu_nominal_height_dist

elseif save_dir == "plots/height_examples/Malibu_fuzzed"
    height_dist = malibu_fuzzed_height_dist
    height_dist_nominal = malibu_nominal_height_dist
end

# Visualize the likelihood of the height distribution
max_height = max(
    mean(height_dist.height_dist) + 3 * std(height_dist.height_dist),
    mean(height_dist_nominal.height_dist) + 3 * std(height_dist_nominal.height_dist)
)

heights = collect(range(0, max_height, length=1000))
likelihoods = height_likelihood(height_dist, heights)
likelihoods_nominal = height_likelihood(height_dist_nominal, heights)

log_likelihoods = height_loglikelihood(height_dist, heights)
log_likelihoods_nominal = height_loglikelihood(height_dist_nominal, heights)

p_l = plot(
    heights, likelihoods, label="Fuzzed", 
    xlabel="Height [m]", ylabel="Likelihood", 
    title="Height Distribution Likelihood", 
    legend=:topright, dpi = 200
)
plot!(
    heights, likelihoods_nominal, label="Nominal"
)
savefig(p_l, joinpath(save_dir, "height_distribution_likelihood.png"))

p_logl = plot(
    heights, log_likelihoods, label="Fuzzed", 
    xlabel="Height [m]", ylabel="Log-Likelihood", 
    title="Height Distribution Log-Likelihood", 
    legend=:topright, dpi = 200
)
plot!(
    heights, log_likelihoods_nominal, label="Nominal"
)
savefig(p_logl, joinpath(save_dir, "height_distribution_loglikelihood.png"))



# Sample from the distribution
n_samples = 1000

# Sample the heights
height_samples = sample_height(height_dist, n_samples)
height_samples_nominal = sample_height(height_dist_nominal, n_samples)

bins = collect(range(0, max_height, length=100))

# Plot the height samples
h = histogram(
    height_samples_nominal, label="Nominal", alpha=0.5,
    bins = bins, normalize=:pdf,
    xlabel="Height [m]", ylabel="Frequency", 
    title="Height Distribution Samples", 
    framestyle = :box, dpi = 200
)
histogram!(
    height_samples, label="Fuzzed", alpha=0.5, 
    bins = bins, normalize=:pdf
)
savefig(h, joinpath(save_dir, "height_distribution_samples.png"))


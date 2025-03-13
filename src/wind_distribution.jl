
using Distributions

# The mixture model type seem so to be 
# MixtureModel{Univariate, Continuous, VonMises{Float64}, Categorical{Float64, Vector{Float64}}}
struct WindDistribution
    wind_dir_dists::MixtureModel #{Distributions.VonMises}
    wind_speed_dist::Normal # To-Do: Change to Gamma Distribution
end

# Constructor
function WindDistribution(
        mixture_dir_means::Vector{Float64}, 
        mixture_dir_κ::Vector{Float64}, 
        mixture_weights::Vector{Float64}, 
        speed_mean::Float64, speed_std::Float64;
        in_deg::Bool=true)

    if in_deg
        mixture_dir_means = [deg2rad(mean) for mean in mixture_dir_means]
    end

    # Wrap the means to be in the range of -π to π
    mixture_means_wrapped = [rem2pi(mean, RoundNearest) for mean in mixture_dir_means]

    mixture_components = [VonMises(mix_mean, mix_κ) 
        for (mix_mean, mix_κ) in zip(mixture_means_wrapped, mixture_dir_κ)]

    # Normalize the mixture_weights
    mixture_weights = mixture_weights ./ sum(mixture_weights)

    return WindDistribution(
        MixtureModel(mixture_components, mixture_weights),
        Normal(speed_mean, speed_std)
    )
end

# Get the standard deviation of each component of the mixture
function dir_stds(wind_dist::WindDistribution, in_deg::Bool=true)
    componend_stds = [std(component) for component in wind_dist.wind_dir_dists.components]
    if in_deg
        return [rad2deg(std(component)) for component in wind_dist.wind_dir_dists.components]
    end
    return componend_stds
end

# Convenience for getting the wind likelihood and loglikelihood that
# wraps around the circle
# The Distributions.jl implementation of the VonMises distribution does not
# wrap around the circle, so we need to get the max of  the likelihoods for the 
# same angle but with branches from -2π to 0 and 0 to 2π
function wind_dir_likelihood_circ(vm::Distributions.VonMises, wind_dir::Float64)
    likelihoods = pdf(vm, [wind_dir - 2π, wind_dir, wind_dir + 2π])
    return maximum(likelihoods)
end

function wind_dir_likelihood_circ(vm::Distributions.VonMises, wind_dirs::Vector{Float64})
    l_below = pdf(vm, wind_dirs .- 2π)
    l_base = pdf(vm, wind_dirs)
    l_above = pdf(vm, wind_dirs .+ 2π)
    
    # Get the element-wise maximum likelihood for each angle
    return [max(l_below[i], l_base[i], l_above[i]) for i in 1:length(wind_dirs)]
end

function wind_dir_likelihood_circ(wind_dist::WindDistribution, wind_dir::Float64)
    likelihoods = pdf(wind_dist.wind_dir_dists, [wind_dir - 2π, wind_dir, wind_dir + 2π])
    return maximum(likelihoods)
end

function wind_dir_likelihood_circ(wind_dist::WindDistribution, wind_dirs::Vector{Float64})
    l_below = pdf(wind_dist.wind_dir_dists, wind_dirs .- 2π)
    l_base = pdf(wind_dist.wind_dir_dists, wind_dirs)
    l_above = pdf(wind_dist.wind_dir_dists, wind_dirs .+ 2π)
    
    return [max(l_below[i], l_base[i], l_above[i]) for i in 1:length(wind_dirs)]
end

# Log likelihoods

function wind_dir_loglikelihood_circ(vm::Distributions.VonMises, wind_dir::Float64)
    loglikelihoods = logpdf(vm, [wind_dir - 2π, wind_dir, wind_dir + 2π])
    return maximum(loglikelihoods)
end

function wind_dir_loglikelihood_circ(vm::Distributions.VonMises, wind_dirs::Vector{Float64})
    logl_below = logpdf(vm, wind_dirs .- 2π)
    logl_base = logpdf(vm, wind_dirs)
    logl_above = logpdf(vm, wind_dirs .+ 2π)

    return [max(logl_below[i], logl_base[i], logl_above[i]) for i in 1:length(wind_dirs)]
end

function wind_dir_loglikelihood_circ(wind_dist::WindDistribution, wind_dir::Float64)
    loglikelihoods = logpdf(wind_dist.wind_dir_dists, [wind_dir - 2π, wind_dir, wind_dir + 2π])
    return maximum(loglikelihoods)
end

function wind_dir_loglikelihood_circ(wind_dist::WindDistribution, wind_dirs::Vector{Float64})
    logl_below = logpdf(wind_dist.wind_dir_dists, wind_dirs .- 2π)
    logl_base = logpdf(wind_dist.wind_dir_dists, wind_dirs)
    logl_above = logpdf(wind_dist.wind_dir_dists, wind_dirs .+ 2π)

    return [max(logl_below[i], logl_base[i], logl_above[i]) for i in 1:length(wind_dirs)]
end

# Convenience for getting the overall likelihood and loglikelihood of the wind
function wind_likelihood(wind_dist::WindDistribution, 
                        wind_dir::Float64, wind_speed::Float64)
    dir_pdf = wind_dir_likelihood_circ(wind_dist, wind_dir)
    speed_pdf = pdf(wind_dist.wind_speed_dist, wind_speed)

    return dir_pdf * speed_pdf
end

function wind_loglikelihood(wind_dist::WindDistribution, 
                        wind_dir::Float64, wind_speed::Float64)
    dir_logpdf = wind_dir_loglikelihood_circ(wind_dist, wind_dir)
    speed_logpdf = logpdf(wind_dist.wind_speed_dist, wind_speed)

    return dir_logpdf + speed_logpdf
end

# Go from wind direction and speed to a wind vector
function to_wind_vec(wind_dir::Float64, wind_speed::Float64, in_deg::Bool=false)
    if in_deg
        wind_dir = deg2rad(wind_dir)
    end
    return wind_speed .* [cos(wind_dir), sin(wind_dir), 0.0]
end

function to_wind_vec(wind_dirs::Vector{Float64}, wind_speeds::Vector{Float64})
    return [to_wind_vec(wd, ws) for (wd, ws) in zip(wind_dirs, wind_speeds)]
end

# Sample from the wind distribution
function sample_wind_dir_speed(wind_dist::WindDistribution, in_deg::Bool=false)
    wind_dir = rand(wind_dist.wind_dir_dists) # Radians
    wind_speed = rand(wind_dist.wind_speed_dist)

    if in_deg
        wind_dir = rad2deg(wind_dir)
    end

    return wind_dir, wind_speed
end

function sample_wind_vec(wind_dist::WindDistribution)
    wind_dir, wind_speed = sample_wind_dir_speed(wind_dist)
    return to_wind_vec(wind_dir, wind_speed)
end

# Sample a trajectory of wind vectors
function sample_wind_vec_trajectory(wind_dist::WindDistribution, n_samples::Int)
    return [sample_wind_vec(wind_dist) for _ in 1:n_samples]
end

function sample_wind_dir_speed_trajectory(wind_dist::WindDistribution, 
        n_samples::Int, in_deg::Bool=false)
    dir_speed_pairs = [sample_wind_dir_speed(wind_dist) for _ in 1:n_samples]

    # Split the pairs
    dirs = [pair[1] for pair in dir_speed_pairs]
    speeds = [pair[2] for pair in dir_speed_pairs]

    # Convert to degrees if needed
    if in_deg
        dirs = rad2deg.(dirs)
    end

    return dirs, speeds
end


# Repeat for a full trajectory of wind speed and direction pairs
function wind_likelihood_traj(wind_dist::WindDistribution, 
                        wind_dirs::Vector{Float64}, wind_speeds::Vector{Float64})
    dir_pdf = wind_dir_likelihood_circ(wind_dist, wind_dirs)
    speed_pdf = pdf(wind_dist.wind_speed_dist, wind_speeds)

    # Return the final product
    return prod(dir_pdf) * prod(speed_pdf)
end

function wind_loglikelihood_traj(wind_dist::WindDistribution, 
                        wind_dirs::Vector{Float64}, wind_speeds::Vector{Float64})
    dir_logpdf = wind_dir_loglikelihood_circ(wind_dist, wind_dirs)
    speed_logpdf = logpdf(wind_dist.wind_speed_dist, wind_speeds)

    # Return the final sum
    return sum(dir_logpdf) + sum(speed_logpdf)
end


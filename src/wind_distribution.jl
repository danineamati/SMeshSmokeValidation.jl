

struct WindDistribution
    wind_dir_dists::MixtureModel{VonMises}
    wind_speed_dist::Normal
end

# Constructor
function WindDistribution(
        mixture_dir_means::Vector{Float64}, 
        mixture_dir_κ::Vector{Float64}, 
        mixture_weights::Vector{Float64}, 
        speed_mean::Float64, speed_std::Float64)

    mixture_components = [VonMises(mixture_dir_means[i], mixture_dir_κ[i]) 
        for i in 1:length(mixture_dir_means)]

    return WindDistribution(
        MixtureModel(mixture_components, mixture_weights),
        Normal(speed_mean, speed_std)
    )
end

# Get the standard deviation of each component of the mixture
function std(wind_dist::WindDistribution)
    return [rad2deg(std(component)) for component in wind_dist.wind_dir_dists.distributions]
end

# Convenience for getting the wind likelihood and loglikelihood that
# wraps around the circle
function wind_dir_likelihood_circ(wind_dist::WindDistribution, wind_dir::Float64)
    # The Distributions.jl implementation of the VonMises distribution does not
    # wrap around the circle, so we need to sum the likelihoods for the same
    # angle but with branches from -2π to 0 and 0 to 2π
    likelihoods = pdf.(wind_dist.wind_dir_dists, [wind_dir, wind_dir + 2π])
    return sum(likelihoods)
end

function wind_dir_loglikelihood_circ(wind_dist::WindDistribution, wind_dir::Float64)
    loglikelihoods = logpdf.(wind_dist.wind_dir_dists, [wind_dir, wind_dir + 2π])
    return sum(loglikelihoods)
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

# Convenience for sampling from the wind distribution
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

function sample_wind_dir_speed_trajectory(wind_dist::WindDistribution, n_samples::Int)
    dir_speed_pairs = [sample_wind_dir_speed(wind_dist) for _ in 1:n_samples]

    # Split the pairs
    dirs = [pair[1] for pair in dir_speed_pairs]
    speeds = [pair[2] for pair in dir_speed_pairs]

    return dirs, speeds
end



# Package all the disturbances into one struct with corresponding functions
# We need to be able to (1) sample a trajectory of the disturbances and 
# (2) calculate the likelihood & log-likelihood of a trajectory of disturbances
# 
# Some parts of the disturbances (i.e., height) only pertains to the start 
# (i.e., initial state distribution) and some parts (i.e., wind) pertains to
# the entire trajectory. 


struct FullDisturbances
    height_dist::HeightDistribution
    wind_dist::WindDistribution
end


struct FullDisturbanceTrajectory
    height::Float64               # Only one height (initial state)
    wind_dirs::Vector{Float64}    # Direction can change
    wind_speeds::Vector{Float64}  # Direction can change
end


# No default constructor. Please use the individual constructors.

# Sample a trajectory of disturbances
function sample_disturbances(disturbances::FullDisturbances, n_timesteps::Int)
    heights = sample_height(disturbances.height_dist, 1)
    wind_dirs, wind_speeds = sample_wind_dir_speed_trajectory(
        disturbances.wind_dist, n_timesteps)

    disturb_traj = FullDisturbanceTrajectory(heights[1], wind_dirs, wind_speeds)
    return disturb_traj
end


# Calculate the likelihood of a trajectory of disturbances
function disturbance_trajectory_likelihood(
        disturbances::FullDisturbances, 
        disturb_traj::FullDisturbanceTrajectory)

    height_pdf = height_likelihood(disturbances.height_dist, 
        disturb_traj.height)
    wind_pdf = wind_likelihood_traj(disturbances.wind_dist, 
        disturb_traj.wind_dirs, disturb_traj.wind_speeds)
    
    return height_pdf * wind_pdf
end


# Calculate the log-likelihood of a trajectory of disturbances
function disturbance_trajectory_log_likelihood(
        disturbances::FullDisturbances, 
        disturb_traj::FullDisturbanceTrajectory)

    height_logpdf = height_loglikelihood(disturbances.height_dist, 
        disturb_traj.height)
    wind_logpdf = wind_loglikelihood_traj(disturbances.wind_dist, 
        disturb_traj.wind_dirs, disturb_traj.wind_speeds)
    
    return height_logpdf + wind_logpdf
end


# Get the most likely disturbance trajectories returned in order of likelihood
function most_likely_disturbance_trajectories(
        loglikelihoods::Vector{Float64}, 
        trajectory_list::Vector{FullDisturbanceTrajectory},
        num_most_likely::Int)

    # Get the most likely trajectories
    most_likely_indices = partialsortperm(
        loglikelihoods, 1:num_most_likely, rev=true)
    most_likely_trajectories = [trajectory_list[i] 
        for i in most_likely_indices]

    return most_likely_trajectories, most_likely_indices
end


function most_likely_disturbance_trajectories(
        disturbances::FullDisturbances, 
        trajectory_list::Vector{FullDisturbanceTrajectory},
        num_most_likely::Int)

    # Log-likelihoods
    loglikelihoods = [
        disturbance_trajectory_log_likelihood(disturbances, traj) 
            for traj in trajectory_list
    ]

    # Get the most likely trajectories and the top indices
    return most_likely_disturbance_trajectories(loglikelihoods, 
                trajectory_list, num_most_likely)
end


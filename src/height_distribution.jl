
using Distributions

# We use Gamma since the height > 0 and should converge to the Normal
# distribution for large Q_release and small h_std
struct HeightDistribution
    height_dist::Gamma
end

# Handle the relationship between Q and h (i.e., units and eventually the 
# Briggs Plume Rise model)
function Q_to_h(Q)
    return Q / 10.0
end

# Constructor
function HeightDistribution(Q_release::Float64, h_std::Float64)
    h_mu = Q_to_h(Q_release)

    gamma_α = h_mu^2 / h_std^2
    gamma_θ = h_std^2 / h_mu
    return HeightDistribution(Gamma(gamma_α, gamma_θ))
end


# Wrap the likelihood and log-likelihood functions in case we change the 
# distribution type in the future
function height_likelihood(height_dist::HeightDistribution, height::Float64)
    return pdf(height_dist.height_dist, height)
end

function height_likelihood(height_dist::HeightDistribution, 
        heights::Vector{Float64})
    return pdf(height_dist.height_dist, heights)
end

# Log-likelihood
function height_loglikelihood(height_dist::HeightDistribution, height::Float64)
    return logpdf(height_dist.height_dist, height)
end

function height_loglikelihood(height_dist::HeightDistribution, 
        heights::Vector{Float64})
    return logpdf(height_dist.height_dist, heights)
end


# Sample from the distribution
function sample_height(height_dist::HeightDistribution, n_samples::Int)
    return rand(height_dist.height_dist, n_samples)
end

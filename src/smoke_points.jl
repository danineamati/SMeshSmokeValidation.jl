# We can sample the smoke from the interior of the burned region at time t.

using LazySets

"""
    sample_smoke_from_poly(burnt_poly::Polygon, n_samples::Int)

Given a polygon that represents the burnt region, we sample n_samples from the
interior of the polygon.
"""
function sample_smoke_from_poly(burnt_poly::Polygon, n_samples::Int;
        k_times::Int=10, verbose::Bool=false)
    # Assert that the polygon is a polygon
    @assert isa(burnt_poly, Polygon)

    # Sample the smoke from the interior of the polygon
    # try k times to sample the smoke, if it fails, return an empty array
    for k_ind in 1:k_times
        try
            # We need to specify LazySets.sample to avoid ambiguity with
            # Distributions and Random
            smoke_samples = LazySets.sample(burnt_poly, n_samples)
            return smoke_samples
        catch e
            if verbose
                println("------> Sampling failed at k_ind = ", k_ind)
            end

            if k_ind == k_times
                println("Error: ", e)
                println("n_samples: ", n_samples)
            end
        end
    end

    return []    
end

"""
    sample_smoke_from_scene(burn_scene_obj::BurnScene, 
                            n_samples::Int, t_ind::Float64)

Given a burn_scene object, we sample the smoke from the interior of the burnt
region at time t.
"""
function sample_smoke_from_scene(burn_scene_obj::BurnScene, 
                                 n_samples::Int, t_ind::Int64;
                                 verbose::Bool=false)
    # If there is no smoke to sample, return an empty array
    if t_ind <= 1
        return []
    end

    smoke_samples = []

    for t_burnt in 1:(t_ind - 1)
        # Sample the smoke from the interior of the polygon at time t_burnt
        # println("Sampling smoke from polygon ", t_burnt)
        push!(smoke_samples, sample_smoke_from_poly(
            burn_scene_obj.burn_polys_t[t_burnt], n_samples))
    end

    # stack them all together as one array
    smoke_samples = vcat(smoke_samples...)

    if verbose
        println("Number of smoke samples: ", length(smoke_samples), 
                    " of ", (t_ind - 1) *n_samples)
    end

    return smoke_samples
end

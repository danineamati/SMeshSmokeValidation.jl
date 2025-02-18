# We can sample the smoke from the interior of the burned region at time t.

using LazySets

"""
    sample_smoke_from_poly(burnt_poly::Polygon, n_samples::Int)

Given a polygon that represents the burnt region, we sample n_samples from the
interior of the polygon.
"""
function sample_smoke_from_poly(burnt_poly::Polygon, n_samples::Int)
    # Sample the smoke from the interior of the polygon
    smoke_samples = sample(burnt_poly, n_samples)
    return smoke_samples
end

"""
    sample_smoke_from_scene(burn_scene_obj::burn_scene, 
                            n_samples::Int, t_ind::Float64)

Given a burn_scene object, we sample the smoke from the interior of the burnt
region at time t.
"""
function sample_smoke_from_scene(burn_scene_obj::burn_scene, 
                                 n_samples::Int, t_ind::Int64)
    # If there is no smoke to sample, return an empty array
    if t_ind <= 1
        return []
    end

    smoke_samples = []

    for t_burnt in 1:(t_ind - 1)
        # Sample the smoke from the interior of the polygon at time t_burnt
        push!(smoke_samples, sample_smoke_from_poly(
            burn_scene_obj.burn_polys_t[t_burnt], n_samples))
    end

    # stack them all together as one array
    smoke_samples = vcat(smoke_samples...)

    return smoke_samples
end

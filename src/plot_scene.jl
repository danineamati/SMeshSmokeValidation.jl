# Given a burn scene and a timestep to plot, we plot the scene at that timestep
# with a consistent burn color scheme as the fire front progresses.

using Plots

"""
    plot_scene(burn_scene_obj::burn_scene, t_ind::Int)

Given a burn scene and a timestep to plot, we plot the scene at that timestep in
red, the scenes that have passed in dark gray, and the scenes that have yet to
come in light gray.
"""
function plot_scene(burn_scene_obj::burn_scene, t_ind::Int;
                    n_smoke_samples::Int=0, 
                    background_image=nothing, 
                    background_x=nothing, background_y=nothing)
    # Extract the polygons and time values from the burn_scene object
    burn_polys_t = burn_scene_obj.burn_polys_t
    t = burn_scene_obj.t

    # Check that the time index is within the bounds of the time vector
    if t_ind > (length(t) + 1) || t_ind < 1
        error("Time index out of bounds")
    end

    # Initialize the plot to the background image if provided
    if !isnothing(background_image)
        p = plot(background_x, background_y,
                 background_image[end:-1:1, :], yflip = false,
                 dpi=300)
    else
        p = plot()
    end

    # Plot the polygons in order
    for ind in eachindex(burn_polys_t)
        if ind < t_ind
            plot!(burn_polys_t[ind], color=:black, alpha=0.5, label="",
                    linealpha=0.0)
        elseif ind == t_ind
            plot!(burn_polys_t[ind], color=:red, label="")
        else
            plot!(burn_polys_t[ind], color=:lightgray, alpha=0.5, label="",
                    linealpha=0.0)
        end
    end

    if n_smoke_samples > 0 && t_ind > 1
        # Sample the smoke from the interior of the polygon at time t
        smoke_samples = sample_smoke_from_scene(burn_scene_obj, n_smoke_samples,
                                                t_ind)
        # The samples are an array of 2-element arrays, so we need to extract
        # the x and y values, which is much easier as a matrix
        smoke_samples_mat = hcat(smoke_samples...)

        # Plot the smoke samples if not empty
        if !isempty(smoke_samples_mat)
            scatter!(smoke_samples_mat[1, :], smoke_samples_mat[2, :],
                color=:lightgray, linecolor=:black, label="Smoke Point Source")
        end
    end

    # Add the axes labels
    xlabel!("Eastings (m)")
    ylabel!("Northings (m)")

    return p

end


"""
    animate_burn_scene(burn_scene_obj::burn_scene)

Given a burn scene, we animate the progression of the fire front as it burns
through the scene using the plot_scene function for each time step.
"""
function animate_burn_scene(burn_scene_obj::burn_scene;
    n_smoke_samples::Int=0, 
    background_image=nothing, 
    background_x=nothing, background_y=nothing)
    # Extract the polygons and time values from the burn_scene object
    burn_polys_t = burn_scene_obj.burn_polys_t

    println("Animating burn scene...")
    println("Number of time steps: ", length(burn_polys_t))
    println("Number of smoke samples: ", n_smoke_samples)
    println("")
    println("burn_polys_t type: ", typeof(burn_polys_t))

    # Plot the polygons in order
    anim = @animate for ind in 1:(1 + length(burn_polys_t))
        plot_scene(burn_scene_obj, ind, 
            n_smoke_samples=n_smoke_samples,
            background_image=background_image, background_x=background_x,
            background_y=background_y)
    end

    return anim
end

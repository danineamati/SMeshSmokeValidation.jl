# Given a burn scene and a timestep to plot, we plot the scene at that timestep
# with a consistent burn color scheme as the fire front progresses.

using Plots

"""
    plot_scene(burn_scene_obj::burn_scene, t_ind::Int)

Given a burn scene and a timestep to plot, we plot the scene at that timestep in
red, the scenes that have passed in dark gray, and the scenes that have yet to
come in light gray.
"""
function plot_scene(burn_scene_obj::burn_scene, t_ind::Int)
    # Extract the polygons and time values from the burn_scene object
    burn_polys_t = burn_scene_obj.burn_polys_t
    t = burn_scene_obj.t

    # Check that the time index is within the bounds of the time vector
    if t_ind > length(t) || t_ind < 1
        error("Time index out of bounds")
    end

    # Initialize the plot
    plot()

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

    # Add the axes labels
    xlabel!("Eastings (m)")
    ylabel!("Northings (m)")

end


"""
    animate_burn_scene(burn_scene_obj::burn_scene)

Given a burn scene, we animate the progression of the fire front as it burns
through the scene using the plot_scene function for each time step.
"""
function animate_burn_scene(burn_scene_obj::burn_scene)
    # Extract the polygons and time values from the burn_scene object
    burn_polys_t = burn_scene_obj.burn_polys_t

    # Plot the polygons in order
    anim = @animate for ind in eachindex(burn_polys_t)
        plot_scene(burn_scene_obj, ind)
    end

    return anim
end

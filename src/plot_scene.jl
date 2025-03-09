# Given a burn scene and a timestep to plot, we plot the scene at that timestep
# with a consistent burn color scheme as the fire front progresses.

using Plots

"""
    plot_scene(burn_scene_obj::BurnScene, t_ind::Int)

Given a burn scene and a timestep to plot, we plot the scene at that timestep in
red, the scenes that have passed in dark gray, and the scenes that have yet to
come in light gray.
"""
function plot_scene(burn_scene_obj::BurnScene, t_ind::Int;
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
                color=:lightgray, linecolor=:black, label="")
                #, label="Smoke Point Source")
        end
    end

    # Add the axes labels
    xlabel!("Eastings (m)")
    ylabel!("Northings (m)")

    return p

end


"""
    animate_burn_scene(burn_scene_obj::BurnScene)

Given a burn scene, we animate the progression of the fire front as it burns
through the scene using the plot_scene function for each time step.
"""
function animate_burn_scene(burn_scene_obj::BurnScene;
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


###############################################
# Plot the smoke plume in global coordinates
###############################################

"""
    plot_single_plume_bounds(plume::PlumeFromPointSource, 
                             bounds::Vector{Float64},
                             wind_vector::Vector{Float64},
                             x_num::Int=101, y_num::Int=103,
                             cmap=:bilbao)

Given a plume object, bounds, and wind vector, we plot the plume in global
coordinates from the Gaussian Plume Model as a density plot over the bounds
provided.

We recommend using relative primes for x_num and y_num to avoid size issues.
"""
function plot_single_plume_bounds(plume::PlumeFromPointSource, 
                                  bounds::Vector{Float64},
                                  wind_vector::Vector{Float64},
                                  x_num::Int=101, y_num::Int=103,
                                  cmap=:bilbao)
    # Extract the bounds
    x_min, x_max, y_min, y_max = bounds

    # For now, we use a single height (future is a DEM)
    height = 0.0

    # Create a meshgrid for the bounds
    x = range(x_min, x_max, length=x_num)
    y = range(y_min, y_max, length=y_num)
    # X = x' .* ones(y_num)
    # Y = ones(x_num)' .* y

    # println("X shape ", size(X))
    # println("Y shape ", size(Y))

    # Compute the plume density over the meshgrid
    density = zeros(x_num, y_num)
    for i in 1:x_num
        for j in 1:y_num
            density[i, j] = query_plume_model(plume, 
                                              [x[i], y[j], height],
                                              wind_vector)
        end
    end

    # Apply a safe log10 to the density
    density[density .<= 1e-10] .= 1e-10
    log10_density = log10.(density)

    # Plot the density
    p = heatmap(x, y, log10_density', aspect_ratio=1, color=cgrad(:bilbao, rev=true),
                xlabel="Eastings (m)", ylabel="Northings (m)")

    # Plot the wind vector
    # quiver!([x_max], [y_max], quiver=(0, 0, wind_vector[1], wind_vector[2]),
    #         color=:black, label="Wind Vector")

    return p
end


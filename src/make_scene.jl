# Given a file that parameterizes the scene and the fire front progression,
# we convert that file to a vertex-based representation of the scene as a 
# function of time.

using LazySets

"""
    burn_scene

A struct that contains the x, y, and t values of the scene.

# Fields
- `burn_polys_t::Vector{Polygon}`: A vector of polygons that represent the
  scene at each time step.
- `t::Vector{Float64}`: A vector of time values that correspond to the scene
  at each time step.
"""
struct burn_scene
    burn_polys_t::Vector{Polygon}
    t::Vector{Float64}
end

"""
    make_burn_polys_from_coords(x::Vector{Float64}, y::Vector{Float64})

Given the x and y coordinates of the scene, we convert the coordinates to a
vertex-based Polygon representation of the scene using the LazySets package.

# Arguments
- `x::Vector{Float64}`: The x-coordinates of the scene.
- `y::Vector{Float64}`: The y-coordinates of the scene.
- `t::Vector{Float64}`: The time values that correspond to the scene at each
  time step.

# Returns
- `burn_polys_t`: A vector of polygons that represent the burn at each time
  step.
"""
function make_burn_polys_from_coords(x, y, t)
    # Initialize the burn_polys_t vector
    burn_polys_t = Vector{Polygon}(undef, length(t))

    # Create the polygons
    for ind in 1:length(t)
        # Extract the x, y values for the polygon
        x_i = x[ind]
        y_i = y[ind]
        # Create the polygon
        poly_i = Polygon([x_i, y_i])
        # Add the polygon to the burn_polys_t vector
        burn_polys_t[ind] = poly_i
    end

    return burn_polys_t
end


"""
    make_scene_from_csv(scene_csv_path::String)

Given a CSV file that parameterizes the scene and the fire front progression,
we convert that file to a vertex-based representation of the scene as a
function of time.

# Arguments
- `scene_csv_path::String`: The path to the CSV file that contains the scene
  parameters.

# Returns
- `burn_scene`: A struct that contains the x, y, and t values of the scene.
"""
function make_scene_from_csv(scene_csv_path::String)
    # Read the CSV file
    df = CSV.read(scene_csv_path)
    # Extract the columns
    x = df[:, :x]
    y = df[:, :y]
    t = df[:, :t]

    burn_polys_t = make_scene_polys_from_coords(x, y, t)

    return burn_scene(burn_polys_t, t)
end

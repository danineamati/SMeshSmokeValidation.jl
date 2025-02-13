# Given a file that parameterizes the scene and the fire front progression,
# we convert that file to a vertex-based representation of the scene as a 
# function of time.

"""
    burn_scene

A struct that contains the x, y, and t values of the scene.

# Fields
- `x::Array{Float64, 1}`: The x values of the scene.
- `y::Array{Float64, 1}`: The y values of the scene.
- `t::Array{Float64, 1}`: The t values of the scene.
"""
struct burn_scene
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    t::Array{Float64, 1}
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

    burn_scene(x, y, t)
end

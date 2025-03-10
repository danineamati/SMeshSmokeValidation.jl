# Given a file that parameterizes the scene and the fire front progression,
# we convert that file to a vertex-based representation of the scene as a 
# function of time.

using LazySets
# using Printf

"""
    BurnScene

A struct that contains the x, y, and t values of the scene.

# Fields
- `burn_polys_t::Vector{Polygon}`: A vector of polygons that represent the
  scene at each time step.
- `t::Vector{Int64}`: A vector of time values that correspond to the scene
  at each time step.
"""
struct BurnScene
    burn_polys_t::Vector{Polygon}
    t::Vector{Int64}
end


"""
    parse_coord_str(coord_str::String)
  
Given a string that contains the coordinates of a polygon, we parse the string
and return a vector of the coordinates. We parse the polygon string at the 
"] [" characters for the coordinate pairs and split the resulting strings at the 
"," characters to get the x and y values.

Example:
```julia
julia> parse_coord_str("[0.0, 0.0] [1.0, 0.0] [1.0, 1.0] [0.0, 1.0]")
4-element Vector{Vector{Float64}}:
 [0.0, 0.0]
 [1.0, 0.0]
 [1.0, 1.0]
 [0.0, 1.0]
```

# Arguments
- `coord_str::String`: A string that contains the coordinates of a polygon.

# Returns
- `coords_poly`: A polygon that represents the coordinates of the input string.
"""
function parse_coord_str_to_poly(coord_str::String)
    return Polygon([
        [parse(Float64, sxy) for sxy in split(s, ',')]
        for s in split(coord_str[2:end-1], "] [")
    ])
end



"""
    make_burn_polys_from_coords(x::Vector{Float64}, y::Vector{Float64})

Given the x and y coordinates of the scene, we convert the coordinates to a
vertex-based Polygon representation of the scene using the LazySets package.

# Arguments
- `x::Vector{Float64}`: The x-coordinates of the scene.
- `y::Vector{Float64}`: The y-coordinates of the scene.
- `t::Vector{Int64}`: The time values that correspond to the scene at each
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


# """
#     make_scene_from_csv(scene_csv_path::String)

# Given a CSV file that parameterizes the scene and the fire front progression,
# we convert that file to a vertex-based representation of the scene as a
# function of time.

# # Arguments
# - `scene_csv_path::String`: The path to the CSV file that contains the scene
#   parameters.

# # Returns
# - `burn_scene`: A struct that contains the x, y, and t values of the scene.
# """
# function make_scene_from_csv(scene_csv_path::String)
#     # Read the CSV file
#     df = CSV.read(scene_csv_path)
#     # Extract the columns
#     x = df[:, :x]
#     y = df[:, :y]
#     t = df[:, :t]

#     burn_polys_t = make_scene_polys_from_coords(x, y, t)

#     return burn_scene(burn_polys_t, t)
# end

"""
Parse a coordinate string of the form "x, y" into a tuple (x, y).

To-Do: Likely we want to return an array rather than a tuple
"""
function parse_coord(str)
    parts = split(str, ",")
    x = parse(Float64, strip(parts[1]))
    y = parse(Float64, strip(parts[2]))
    return (x, y)
end

# r""" """ is necessary due to the escape characters in the regex
r"""
    parse_line_to_poly(coord_regex=r"\[([\d\.\-+eE,\s]+)\]")

Given a line from a file that contains the coordinates of a polygon, we parse 
the line and return a polygon.

The coordinate regex is set to r"\[([\d\.\-+eE,\s]+)\]" by default, which
matches the "text [x1, y1] [x2, y2] ..." format.

To-Do: Align better with the parse_coord_str_to_poly function.
"""
function parse_line_to_poly(lines, coord_regex=r"\[([\d\.\-+eE,\s]+)\]")
  polys_per_line = Polygon[]
  for line in lines[3:end]
      matches = collect(eachmatch(coord_regex, line))
      coords = [collect(parse_coord(m.captures[1])) for m in matches]
      poly = Polygon(coords)
      push!(polys_per_line, poly)
  end
  return polys_per_line
end


function load_burn_area_file(filename)
    lines = readlines(filename)
    # Parse the coordinates into polygons
    polys_per_line = parse_line_to_poly(lines)

    # Get the JSON name (i.e., area3.json indicates time 3)
    # The first row is a header and the second row is the total area
    burn_area_nums = [parse(Int, filter(isdigit, collect(split(l)[1]))[1]) 
                        for l in lines[3:end]]
    
    # Reorder with the area numbers
    polys_per_line = [polys_per_line[i] for i in burn_area_nums]

    # Convert to a BurnScene object
    burn_scene_obj = BurnScene(polys_per_line, 1:length(polys_per_line))
    return burn_scene_obj
end




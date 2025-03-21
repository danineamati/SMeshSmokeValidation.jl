# Given a file that parameterizes the scene and the fire front progression,
# we convert that file to a vertex-based representation of the scene as a 
# function of time.

using LazySets
using Printf
using Plots  

const COORD_REGEX = r"\[([\d\.\-+eE,\s]+)\]"

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
    snode_locations::Vector{Vector{Float64}}
    reference_bounds_x::Vector{Float64}
    reference_bounds_y::Vector{Float64}
    perimeter_sample_points::Vector{Vector{Float64}}
end

# Default constructor with just the burn polygons and time values
BurnScene(burn_polys_t::Vector{Polygon}, t::Vector{Int64}) = BurnScene(
    burn_polys_t, t, 
    Vector{Vector{Float64}}(), Vector{Float64}(), Vector{Float64}(),
    Vector{Vector{Float64}}()
)

# -----------------------------
# Parsing and distance functions
# -----------------------------
"""
Parse a coordinate string of the form "x, y" into a tuple (x, y).
"""
function parse_coord(str)
    parts = split(str, ",")
    x = parse(Float64, strip(parts[1]))
    y = parse(Float64, strip(parts[2]))
    return (x, y)
end

dist(p1, p2) = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)

# -----------------------------
# Load coordinates and get lengths
# -----------------------------

"""
Reads the file `filename`, extracts all coordinate pairs, and returns:
  - The closed coordinate list (with coords[end] = coords[1] if close_polygon=true)
  - A vector of segment lengths between consecutive coords.
"""
function get_coords_and_lengths(filename; close_polygon::Bool=true)
    txt = read(filename, String)
    matches = collect(eachmatch(COORD_REGEX, txt))
    coords = [parse_coord(m.captures[1]) for m in matches]

    # Close the polygon if not already closed
    if close_polygon && !isempty(coords) && coords[1] != coords[end]
        push!(coords, coords[1])
    end

    lengths = [dist(coords[i], coords[i+1]) for i in 1:length(coords)-1]
    return coords, lengths
end

"""
Build a cumulative distance array where cumdist[i] = total distance
from the first coordinate up to coords[i].
  - cumdist[1] = 0.0
  - cumdist[end] = total perimeter
"""

function build_cumulative_distances(coords, lengths)
    cumdist = Vector{Float64}(undef, length(coords))
    cumdist[1] = 0.0
    for i in 2:length(coords)
        cumdist[i] = cumdist[i-1] + lengths[i-1]
    end
    return cumdist
end

# -----------------------------
# Sampling function
# -----------------------------
"""
Sample n points (by default n=100) uniformly around the perimeter.
  - coords: closed coordinate list
  - lengths: vector of segment lengths
  - n: number of samples

Returns an array of sampled point tuples [(x1, y1), (x2, y2), ...].
"""
function sample_perimeter_points(coords, lengths; n=100)
    cumdist = build_cumulative_distances(coords, lengths)
    total_perim = last(cumdist)

    # Distances at which we sample (0 up to just before total_perim, in n equal increments)
    sample_dists = [i*(total_perim/n) for i in 0:(n-1)]

    sampled_points = Vector{Tuple{Float64, Float64}}(undef, n)
    seg_idx = 2  # segment index for coords/cumdist

    for (j, d) in enumerate(sample_dists)
        # Move seg_idx forward until d lies in [cumdist[seg_idx-1], cumdist[seg_idx])
        while seg_idx <= length(cumdist) && d > cumdist[seg_idx]
            seg_idx += 1
        end
        if seg_idx > length(cumdist)
            seg_idx = length(cumdist)
        end

        # Now interpolate between coords[seg_idx-1] and coords[seg_idx]
        d0 = cumdist[seg_idx-1]
        d1 = cumdist[seg_idx]
        t = (d - d0)/(d1 - d0)  # fraction along that segment

        p1 = coords[seg_idx-1]
        p2 = coords[seg_idx]

        x = p1[1] + t*(p2[1] - p1[1])
        y = p1[2] + t*(p2[2] - p1[2])
        sampled_points[j] = (x, y)
    end

    return sampled_points
end

# -----------------------------
# New Function: Sample Perimeter from File
# -----------------------------
"""
sample_perimeter_from_file(filename; n=100)

Takes a file containing total perimeter positions (coordinate pairs) and returns
a Vector{Vector{Float64}}, where each inner vector represents a sampled point [x, y].
"""
function sample_perimeter_from_file(filename::String; n::Int=100)::Vector{Vector{Float64}}
    coords, lengths = get_coords_and_lengths(filename)
    sampled_tuples = sample_perimeter_points(coords, lengths; n=n)
    # Convert each tuple (x, y) to a Vector{Float64}
    perimeter_sample_points = [ [pt[1], pt[2]] for pt in sampled_tuples ]
    return perimeter_sample_points
end

"""
[DEPRECATED]

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
[DEPRECATED]

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

# r""" """ is necessary due to the escape characters in the regex
r"""
    parse_lines_to_coords(coord_regex=r"\[([\d\.\-+eE,\s]+)\]")
    parse_lines_to_coords(coord_regex=r"\[([\d\.\-+eE,\s]+)\]")

Given a line from a file that contains the coordinates of a polygon, we parse 
the line and return a vector of coordinates.
the line and return a vector of coordinates.

The coordinate regex is set to r"\[([\d\.\-+eE,\s]+)\]" by default, which
matches the "text [x1, y1] [x2, y2] ..." format.

To-Do: Align better with the parse_coord_str_to_poly function.
"""
function parse_lines_to_coords(lines; coord_regex=r"\[([\d\.\-+eE,\s]+)\]",
    start_line=3)  
  coords_per_line = []
  for line in lines[start_line:end]
      matches = collect(eachmatch(coord_regex, line))
      coords = [collect(parse_coord(m.captures[1])) for m in matches]
      push!(coords_per_line, coords)
  end
  return coords_per_line
end


function parse_lines_to_poly(lines, coord_regex=r"\[([\d\.\-+eE,\s]+)\]")
  polys_per_line = [Polygon(coords) 
      for coords in parse_lines_to_coords(lines, coord_regex=coord_regex)]
  return polys_per_line
end


function load_burn_area_file(filename)
  lines = readlines(filename)
  # Parse the coordinates into polygons
  polys_per_line = parse_lines_to_poly(lines)
  lines = readlines(filename)
  # Parse the coordinates into polygons
  polys_per_line = parse_lines_to_poly(lines)

  # Get the JSON name (i.e., area3.json indicates time 3)
  # The first row is a header and the second row is the total area
  burn_area_nums = [parse(Int, filter(isdigit, collect(split(l)[1]))[1]) 
                      for l in lines[3:end]]
  
  # Reorder with the area numbers
  polys_per_line = [polys_per_line[i] for i in burn_area_nums]

  return polys_per_line
end


function load_locations_file(filename; start_line=2)
  lines = readlines(filename)
  # Parse the coordinates
  locations = parse_lines_to_coords(lines, start_line=start_line)
  # Currently, there is an extract set of brackets around the snode locations
  # (each location is one point, not a polygon)
  locations = [s[1] for s in locations]

  return locations
end


function load_burn_scene_from_files(
      burn_area_filename, 
      snode_locations_filename,
      reference_bounds_filename = nothing,
      total_burn_area_filename = nothing)
  burn_polys = load_burn_area_file(burn_area_filename)
  snode_locations = load_locations_file(snode_locations_filename)

  if isnothing(total_burn_area_filename)
    sample_locations = []
  else
    sample_locations = sample_perimeter_from_file(total_burn_area_filename)
  end

  if isnothing(reference_bounds_filename)
    reference_bounds_x = []
    reference_bounds_y = []
  else
    reference_bounds = load_locations_file(reference_bounds_filename)

    # Sort to min, max
    reference_bounds_x = sort([r[1] for r in reference_bounds])
    reference_bounds_y = sort([r[2] for r in reference_bounds])

  end

  burn_scene_obj = BurnScene(burn_polys, 1:length(burn_polys), 
      snode_locations, reference_bounds_x, reference_bounds_y, sample_locations)
  return burn_scene_obj
end


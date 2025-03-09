#!/usr/bin/env julia

using Printf
using Plots  # For plotting

const COORD_REGEX = r"\[([\d\.\-+eE,\s]+)\]"

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

    # Distances at which we sample
    # (0 up to just before total_perim, in n equal increments)
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

        # Now we interpolate between coords[seg_idx-1] and coords[seg_idx]
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
# Plotting
# -----------------------------
"""
Plot the polygon and the sampled points.
"""
function plot_polygon_and_samples(coords, sampled_points)
    # Separate x and y for original polygon
    xcoords = [p[1] for p in coords]
    ycoords = [p[2] for p in coords]

    # x/y for sampled points
    xsamples = [p[1] for p in sampled_points]
    ysamples = [p[2] for p in sampled_points]

    # Plot the polygon as a path
    plt = plot(xcoords, ycoords, seriestype=:path, label="Polygon", aspect_ratio=:equal)

    # Plot the sampled points on top
    scatter!(plt, xsamples, ysamples, label="Samples")

    return plt
end

# -----------------------------
# Main driver
# -----------------------------
function main()
    if length(ARGS) < 1
        println("Usage: julia Sample_Perimeter.jl /path/to/TotalArea.txt [n=100]")
        return
    end

    filename = ARGS[1]
    n = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100

    coords, lengths = get_coords_and_lengths(filename)
    perimeter = sum(lengths)
    @printf("Total perimeter = %.6f\n", perimeter)

    # Print segment lengths
    for (i, l) in enumerate(lengths)
        @printf("Segment %d length = %.6f\n", i, l)
    end

    # Sample n points
    sampled_points = sample_perimeter_points(coords, lengths; n=n)
    println("\nSampled points (n = $n):")
    for (i, pt) in enumerate(sampled_points)
        @printf("  %3d: (%.6f, %.6f)\n", i, pt[1], pt[2])
    end

    # Plot it
    plt = plot_polygon_and_samples(coords, sampled_points)
    display(plt)  # Show in a window or notebook (depends on environment)

    # Uncomment below to save the plot to a file, e.g. PNG:
    # savefig(plt, "polygon_samples.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

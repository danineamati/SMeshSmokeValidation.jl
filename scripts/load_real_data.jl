
using LazySets
using Images
using TiffImages
using Random

using SMeshSmokeValidation

# Include the functions defined in 'sample_perimeter.jl' here.
include("sample_perimeter.jl")

# Set the RNG seed for reproducibility
Random.seed!(42)

dataset = "HenryCoe"

# Save directory for the plots
save_dir = "plots/" * dataset
if !isdir(save_dir)
    mkdir(save_dir)
end

# Pull the background image

background_img = TiffImages.load("data\\BurnData\\DEMs_and_Buffered_Burns\\DEM_" * dataset * ".tif")
println("Background image size: ", size(background_img))
println("Background image type: ", typeof(background_img))
println("Background image color type: ", eltype(background_img))

println(background_img[10, 3])

# Convert type from Q0f15 to Float32 (i.e., 255 levels of gray)
background_img = Gray{Float32}.(background_img)
#reinterpret(Float32, background_img)
#convert(Array{Float32}, background_img)
println("Background image size: ", size(background_img))
println("Background image type: ", typeof(background_img))
println("Background image color type: ", eltype(background_img))

println(background_img[10, 3])

println(background_img[end:-1:1, :][10, 3])



# error("Stop here")

references = "data\\BurnData\\" * dataset * "\\ReferencePoints\\referencelocations.txt"
references_lines = readlines(references)

ref_points = Vector{Float64}[]
for ref_line in references_lines[2:end]
    matches = collect(eachmatch(COORD_REGEX, ref_line))
    coords = [collect(parse_coord(m.captures[1])) for m in matches]
    push!(ref_points, coords[1])
end

background_img_x_bound = [ref_points[1][1], ref_points[2][1]]
background_img_y_bound = [ref_points[1][2], ref_points[2][2]]

p_back = plot(background_img_x_bound, background_img_y_bound,
                 background_img[end:-1:1, :], yflip = false,
                 dpi=300)
savefig(p_back, save_dir * "/background_image_test.png")

# println("Background image bounds: ", background_img_x_bound, background_img_y_bound)

# Parse the burn areas

burnareas_filename = "data\\BurnData\\" * dataset * "\\BurnAreas\\" * dataset * "BurnAreas.txt"

lines = readlines(burnareas_filename)

# Get the JSON name (first text before the spaces)
json_name = split(lines[end])[1]

burn_area_nums = [parse(Int, filter(isdigit, collect(split(l)[1]))[1]) for l in lines[3:end]]
# println("Area numbers: ", burn_area_nums)
num_steps = length(burn_area_nums)


# Get the coordinates
polys_per_line = Polygon[]
for line in lines[3:end]
    matches = collect(eachmatch(COORD_REGEX, line))
    coords = [collect(parse_coord(m.captures[1])) for m in matches]
    poly = Polygon(coords)
    push!(polys_per_line, poly)
end

# Reorder with the area numbers
# println("Burn Area numbers: ", burn_area_nums)
polys_per_line = [polys_per_line[i] for i in burn_area_nums]


# Convert to a BurnScene object
burn_scene_obj = BurnScene(polys_per_line, 1:num_steps)


for t_ind in 1:(num_steps+1)
    # println("Plotting burn scene at t = ", t_ind)
    p = plot_scene(burn_scene_obj, t_ind, 
        n_smoke_samples=10)
    
    # Save the plot to file 
    savefig(p, save_dir * "/burn_scene_henry_coe_t$(t_ind).png")

end


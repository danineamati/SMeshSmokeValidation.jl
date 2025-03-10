
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

reference_bounds_filename = "data\\BurnData\\" * dataset * "\\ReferencePoints\\referencelocations.txt"
references_lines = readlines(reference_bounds_filename)

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
snode_locations_filename = "data\\BurnData\\" * dataset * "\\SNodeLocations\\snodelocations.txt"

burn_scene_obj = load_burn_scene_from_files(
    burnareas_filename, snode_locations_filename, reference_bounds_filename)

println("Burn Scene Specifications:")
println("Number of burn polygons: ", length(burn_scene_obj.burn_polys_t))
println("Number of time steps: ", length(burn_scene_obj.t))
println("Number of snode locations: ", length(burn_scene_obj.snode_locations))
println("Reference Bounds in x: ", burn_scene_obj.reference_bounds_x)
println("Reference Bounds in y: ", burn_scene_obj.reference_bounds_y)

for t_ind in burn_scene_obj.t
    # println("Plotting burn scene at t = ", t_ind)
    p = plot_scene(burn_scene_obj, t_ind, 
        n_smoke_samples=10)
    
    # Save the plot to file 
    savefig(p, save_dir * "/burn_scene_henry_coe_t$(t_ind).png")

end


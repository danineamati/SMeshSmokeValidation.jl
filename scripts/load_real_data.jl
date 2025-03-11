using LazySets
using Images
using FileIO
using Random
using SMeshSmokeValidation
using Colors
using CSV, DataFrames

include("sample_perimeter.jl")
Random.seed!(42)

# Model Input Paramters
dataset = "HenryCoe"      # Desired burn (e.g., "HenryCoe", "Shasta", "Malibu")
simulation = 1            # 1-3
moisture_level = "moist"  # "moderate", "moist", "wet"

# Save directory for the plots
save_dir = joinpath("plots", dataset)
if !isdir(save_dir)
    mkdir(save_dir)
end

# Pull the background image
img_path = joinpath("data", "BurnData", "DEMs_and_Buffered_Burns", "Background_$(dataset).png")
background_img = load(img_path)
background_img = RGB.(background_img)

println("Background image size: ", size(background_img))
println("Background image type: ", typeof(background_img))
println("Background image color type: ", eltype(background_img))
println("A sample pixel: ", background_img[10, 3])

reference_bounds_filename = joinpath("data", "BurnData", dataset, "ReferencePoints", "referencelocations.txt")
references_lines = readlines(reference_bounds_filename)

# burnareas_filename
burnareas_filename = joinpath("data", "BurnData", dataset, "BurnAreas", "$(dataset)BurnAreas.txt")
snode_locations_filename = joinpath("data", "BurnData", dataset, "SNodeLocations", "snodelocations_sim$(simulation).txt")

burn_scene_obj = load_burn_scene_from_files(
    burnareas_filename, snode_locations_filename, reference_bounds_filename)

println("Burn Scene Specifications:")
println("Number of burn polygons: ", length(burn_scene_obj.burn_polys_t))
println("Number of time steps: ", length(burn_scene_obj.t))
println("Number of snode locations: ", length(burn_scene_obj.snode_locations))
println("Reference Bounds in x: ", burn_scene_obj.reference_bounds_x)
println("Reference Bounds in y: ", burn_scene_obj.reference_bounds_y)

# Q values for the burn scene:
Acreage_Consumption = joinpath("data", "FuelConsumptionProfiles", "$(dataset)_Scenarios.csv")
df = CSV.read(Acreage_Consumption, DataFrame)

# Draw key environmental values
consumed_total = parse(Float64, df[11, Symbol(moisture_level)])  

duration = 3600*24                # seconds of one day burn
Q = consumed_total*1e6 / duration # grams / second
println("Q value (g/acre/sec): ", Q)

# Extract coordinate extents from your burn scene reference bounds
# (Assuming burn_scene_obj is already loaded)
xlims = burn_scene_obj.reference_bounds_x  # e.g., [342384.49223234225, 343981.7686972534]
ylims = burn_scene_obj.reference_bounds_y  # e.g., [3.768196725691409e6, 3.769207180586237e6]

# Create coordinate arrays that span the image dimensions
nx = size(background_img, 2)  # number of columns (x-dimension)
ny = size(background_img, 1)  # number of rows (y-dimension)
background_x = range(xlims[1], stop=xlims[2], length=nx)
background_y = range(ylims[1], stop=ylims[2], length=ny)

num_steps = length(burn_scene_obj.t) + 1

# Now call plot_scene with the background image and its coordinates
for t_ind in 1:(length(burn_scene_obj.t)+1)
    p = plot_scene(burn_scene_obj, t_ind, 
                   n_smoke_samples=10, 
                   background_image=background_img,
                   background_x=background_x,
                   background_y=background_y)
    
    # Save the plot to file 
    savefig(p, joinpath(save_dir, "burn_scene_test_$(dataset)_topobkgd_$(t_ind).png"))
end

println("Done!")

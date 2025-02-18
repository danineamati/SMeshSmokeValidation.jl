using LazySets
using Plots

using SMeshSmokeValidation

# Test that we can call objects from SMeshSmokeValidation
println("SMeshSmokeValidation loaded successfully")

P1 = Polygon([[0.0, 0], [0, 2], [2, 2], [2, 0], [1, 1]]);
P2 = Polygon([[0.0, -1], [0, 0], [1, 1], [2, 0], [3, -1], [1, -1]]);
P3 = Polygon([[0.0, -2], [0, -1], [3, -1], [3, -2], [2, -3]]);
P4 = Polygon([[0.0, -3], [0, -2], [2, -3], [3, -4], [1, -4]]);
P_list = [P1, P2, P3, P4]

burn_scene_obj = burn_scene(P_list, Vector{Float64}(1:length(P_list)))

# Plot the burn scene
plot_scene(burn_scene_obj, 3, n_smoke_samples=10)

# Animate the burn scene
anim = animate_burn_scene(burn_scene_obj, n_smoke_samples=100)
gif(anim, "plots/burn_progression_example.gif", fps = 1)

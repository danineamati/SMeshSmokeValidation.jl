# We use a Gaussian plume model to simulate the smoke dispersion. The model is
# based on the following assumptions:
#
# 1. The smoke is emitted from a point source.
# 2. The smoke is emitted at a constant rate.
# 3. The smoke is emitted at a constant wind speed and direction.

using LinearAlgebra

# Vertical column parameters (see Wikipedia)
# https://en.wikipedia.org/wiki/Atmospheric_dispersion_modeling
# Using Julia matrix notation
VERTICAL_COLUMN_VALUES_TABLE = [
    -1.104  0.9878  -0.0076   4.679  -1.7172   0.2770;
    -1.634  1.0350  -0.0096  -1.999   0.8752   0.0136;
    -2.054  1.0231  -0.0076  -2.341   0.9477  -0.0020;
    -2.555  1.0423  -0.0087  -3.186   1.1737  -0.0316;
    -2.754  1.0106  -0.0064  -3.783   1.3010  -0.0450;
    -3.143  1.0148  -0.0070  -4.490   1.4024  -0.0540
]

"""
    atmo_col_params

The parameters (i,j,k) for y and z for the air column stability classes.
"y" is the cross-wind direction and "z" is the vertical direction.
"""
struct atmo_col_params
    iy::Float64
    jy::Float64
    ky::Float64
    iz::Float64
    jz::Float64
    kz::Float64
end

# Enumerate the air column stability classes as a dictionary
AIR_COLUMN_STABILITY_CLASSES = Dict(
    "A" => atmo_col_params(VERTICAL_COLUMN_VALUES_TABLE[1, :]...),
    "B" => atmo_col_params(VERTICAL_COLUMN_VALUES_TABLE[2, :]...),
    "C" => atmo_col_params(VERTICAL_COLUMN_VALUES_TABLE[3, :]...),
    "D" => atmo_col_params(VERTICAL_COLUMN_VALUES_TABLE[4, :]...),
    "E" => atmo_col_params(VERTICAL_COLUMN_VALUES_TABLE[5, :]...),
    "F" => atmo_col_params(VERTICAL_COLUMN_VALUES_TABLE[6, :]...)
)

"""
    σChange(x_local::Float64, iyz::Float64, jyz::Float64, kyz::Float64)

The change in the standard deviation of the Gaussian plume model as the smoke
plume disperses downwind.

The parameters iyz, jyz, and kyz are the coefficients of the polynomial that
come from the air column stability class.
"""
function σChange(x_local::Float64, iyz::Float64, jyz::Float64, kyz::Float64)
    return exp(iyz + jyz * log(x_local) + kyz * log(x_local)^2)
end


"""
    cross_wind_dispersion(x_local::Float64, y_local::Float64, σ_y::Float64)

The cross-wind dispersion of the Gaussian plume model.
"""
function cross_wind_dispersion(x_local::Float64, y_local::Float64, σ_y::Function)
    σy_curr = σ_y(x_local)
    return exp(- y_local^2 / (2 * σy_curr^2)) / (σy_curr * sqrt(2 * π))
end


"""
    along_wind_dispersion(x_local::Float64, z_local::Float64, σ_z::Float64,
                          h_center::Float64)

The along-wind dispersion of the Gaussian plume model.
"""
function along_wind_dispersion(x_local::Float64, z_local::Float64, σ_z::Function,
                               h_center::Float64)
    σz_curr = σ_z(x_local)
    return exp(- (z_local - h_center)^2 / (2 * σz_curr^2)) / (σz_curr * sqrt(2 * π))
end


"""
    plume_model_point_estimate(
                x_local::Float64, y_local::Float64, z_local::Float64, 
                Q::Float64, u_wind_local::Float64, h_local::Float64,
                σ_y::Function, σ_z::Function)

The point estimate of the Gaussian plume model at a point (x_local, y_local,
z_local) in the smoke plume.
"""
function plume_model_point_estimate(
    x_local::Float64, y_local::Float64, z_local::Float64, 
    Q::Float64, u_wind_local::Float64, h_local::Float64,
    σ_y::Function, σ_z::Function)

    # if the query is upwind of the source, return 0
    if x_local < 0
        return 0
    end

    cross_wind_curr = cross_wind_dispersion(x_local, y_local, σ_y)
    along_wind_curr = along_wind_dispersion(x_local, z_local, σ_z, h_local)

    point_estimate = (Q / u_wind_local) * cross_wind_curr * along_wind_curr

    # println("Point estimate: ", point_estimate)

    return point_estimate
end



#############################
# Package into a struct
#############################

"""
    PlumeFromPointSource

A struct that represents a Gaussian plume model from a point source, enabling
query without having to repass the parameters.
"""
struct PlumeFromPointSource
    Q::Float64
    h::Float64
    σ_y::Function
    σ_z::Function
    source_xyz::Vector{Float64}
end


"""
    PlumeFromPointSource(Q::Float64, h::Float64, stability_class::String
                         source_xyz::Vector{Float64})

Construct a PlumeFromPointSource object.
"""
function PlumeFromPointSource(Q::Float64, h::Float64, stability_class::String,
                              source_xyz::Vector{Float64})
    
    atmo_params = AIR_COLUMN_STABILITY_CLASSES[stability_class]
    σ_y(x_local) = σChange(x_local, atmo_params.iy, atmo_params.jy, atmo_params.ky)
    σ_z(x_local) = σChange(x_local, atmo_params.iz, atmo_params.jz, atmo_params.kz)

    return PlumeFromPointSource(Q, h, σ_y, σ_z, source_xyz)
end


"""
    global_to_local_coords(xyz_local::Vector{Float64}, 
                           xyz_global::Vector{Float64},
                           wind_global::Vector{Float64})

Get a query point and wind vector in global coordinates and return them in
local coordinates.

The wind in (u, v, w) will become (1, 0, 0) in local directions (i.e., along x).
"""
function global_to_local_coords(xyz_local::Vector{Float64}, 
                                xyz_global::Vector{Float64},
                                wind_global::Vector{Float64})
    # Shift the global coordinates to the local coordinates
    xyz_centered = xyz_global - xyz_local

    # Get the local coordinate system
    wind_speed = norm(wind_global)
    # This is the x direction in the local coordinates
    wind_x_local = wind_global / wind_speed 
    # z direction is the same as the global z
    z_local = [0, 0, 1]
    # y direction is the cross product of z and x
    y_local = cross(z_local, wind_x_local)
    # Normalize the y direction (should be normalized, but just in case)
    y_local = y_local / norm(y_local)

    # Get the local coordinates
    xyz_local = [dot(xyz_centered, wind_x_local), 
                 dot(xyz_centered, y_local), 
                 dot(xyz_centered, z_local)]

    # Check that the magnitude was preserved
    @assert isapprox(norm(xyz_centered), norm(xyz_local))
    
    return xyz_local
end


"""
    global_to_local_coords(plume::PlumeFromPointSource, 
                            xyz_global::Vector{Float64},
                            wind_global::Vector{Float64})

Get a query point and wind vector in global coordinates and return them in
local coordinates.
"""
function global_to_local_coords(plume::PlumeFromPointSource, 
                                xyz_global::Vector{Float64},
                                wind_global::Vector{Float64})
    return global_to_local_coords(plume.source_xyz, xyz_global, wind_global)    
end



"""
    query_plume_model(plume_model::PlumeFromPointSource, 
                      xyz_global::Vector{Float64}, wind_global::Vector{Float64})

Query the Gaussian plume model at a point in global coordinates knowing the 
wind in global coordinates.
"""
function query_plume_model(plume_model::PlumeFromPointSource, 
                           xyz_global::Vector{Float64}, 
                           wind_global::Vector{Float64})
    # Get the local coordinates
    xyz_local = global_to_local_coords(plume_model, xyz_global, wind_global)

    # Get the point estimate
    return plume_model_point_estimate(xyz_local..., 
                    plume_model.Q, norm(wind_global), plume_model.h,
                    plume_model.σ_y, plume_model.σ_z)
end



#############################
# Querying multiple plumes
#############################


"""
    query_multiple_plumes_points(
            plumes::Vector{PlumeFromPointSource}, 
            xyz_global_points::Vector{Vector{Float64}}, 
            wind_global::Vector{Float64})

Query the Gaussian plume model at multiple points in global coordinates knowing 
the wind in global coordinates.

The code is not particularly vectorized, so it is relying on Julia's JIT for
performance.
"""
function query_multiple_plumes_points(
        plumes::Vector{PlumeFromPointSource}, 
        xyz_global_points::Vector{Vector{Float64}}, 
        wind_global::Vector{Float64})
    
    # Initialize the density
    density_vec = zeros(length(xyz_global_points))

    # Query the plume model at each point
    for (pt_ind, xyz_global) in enumerate(xyz_global_points)
        for plume in plumes
            density_vec[pt_ind] += query_plume_model(
                plume, xyz_global, wind_global)
        end
    end

    return density_vec
end


using DrWatson
@quickactivate "GyrotronMagnet"

using GyrotronMagnet
using Gridap
using Plots

"""
Example script demonstrating the magnetostatic FEM workflow.

This script shows how to:
1. Define coil geometry
2. Load or generate mesh
3. Create current density field
4. Solve magnetostatic problem
5. Extract and visualize magnetic field
"""

println("=" ^ 60)
println("GyrotronMagnet - Magnetostatic Example")
println("=" ^ 60)

# Step 1: Define coil parameters
println("\n1. Defining coil geometry...")
# create_coil_geometry(r, z, I, n_turns; r_inner, r_outer, z_length)
# Note: r_inner and r_outer are ABSOLUTE radii (not relative to r)
coil1 = create_coil_geometry(0.1, 0.5, 1000.0, 100;
    r_inner = 0.08,    # Inner radius [m]
    r_outer = 0.12,    # Outer radius [m]
    z_length = 0.1     # Axial length [m]
)

coil2 = create_coil_geometry(0.1, 0.7, 1000.0, 100;
    r_inner = 0.08,
    r_outer = 0.12,
    z_length = 0.1
)

coils = [coil1, coil2]
println("   Created $(length(coils)) coils")

# Step 2: Create mesh from coil geometry
println("\n2. Creating mesh from coil geometry...")
domain = (0.0, 0.3, 0.0, 1.0)  # (r_min, r_max, z_min, z_max)
resolution = (120, 180)
model = CartesianDiscreteModel(domain, resolution)
labels = get_face_labeling(model)
println("   Mesh generated on domain r∈[0,0.3], z∈[0,1.0]")

# Step 3: Create current density field (scalar J_φ for axisymmetric)
println("\n3. Creating current density field...")
J_φ = create_current_density(coils, model)
println("   Current density field created (azimuthal component)")

# Step 4: Solve magnetostatic problem (2D axisymmetric)
println("\n4. Solving magnetostatic problem (2D axisymmetric)...")
println("   This may take a few moments...")
# Note: Boundary tag handling - try to add "outer_boundary" tag if using Gridap functions
# For Gmsh meshes, tags should already be defined in the mesh file
try
    using Gridap: add_tag_from_tags!
    # Try to identify boundary tags (common for rectangular domains: 1,2,3,4 or 6,7,8,9)
    boundary_tags = [1, 2, 3, 4]
    add_tag_from_tags!(labels, "outer_boundary", boundary_tags)
catch
    # Tag might already exist or mesh has different tag structure
    # Will use default boundary handling in solver
end
solution = solve_magnetostatic(model, J_φ)  # Uses default μ₀ = 4π×10⁻⁷
println("   Solution completed")

# Step 5: Extract field on axis
println("\n5. Extracting magnetic field on axis...")
z_points = collect(0.0:0.01:1.0)
B_z_values = extract_field_on_axis(solution.B_z, z_points)
println("   Field extracted at $(length(z_points)) points")

# Step 6: Visualize results
println("\n6. Visualizing results...")
p = plot(z_points, B_z_values, 
         xlabel="z [m]", 
         ylabel="B_z [T]", 
         title="Magnetic Field on Axis (2D Axisymmetric)",
         linewidth=2,
         legend=false)
display(p)

# Save plot
plot_path = plotsdir("field_on_axis.png")
mkpath(dirname(plot_path))
savefig(p, plot_path)
println("   Plot saved to: $plot_path")

# Step 7: Save results
println("\n7. Saving results...")
params = Dict(
    "coils" => length(coils),
    "coil1_r" => coil1.r,
    "coil1_z" => coil1.z,
    "coil1_I" => coil1.I,
    "coil2_r" => coil2.r,
    "coil2_z" => coil2.z,
    "coil2_I" => coil2.I,
)

save_results(solution.A_φ, solution.B_z, params, "magnetostatic_solution")
println("   Results saved")

println("\n" * "=" ^ 60)
println("Example completed successfully!")
println("=" ^ 60)


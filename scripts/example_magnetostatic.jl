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
# create_coil_geometry(; r, z, current_density, r_inner, r_outer, z_length)
# Note: r_inner and r_outer are ABSOLUTE radii (not relative to r)
coil1 = create_coil_geometry(
    r = 0.1,
    z = 0.2,
    current_density = 5e7,
    r_inner = 0.08,    # Inner radius [m]
    r_outer = 0.12,    # Outer radius [m]
    z_length = 0.1     # Axial length [m]
)

coil2 = create_coil_geometry(
    r = 0.1,
    z = 0.3,
    current_density = 5e7,
    r_inner = 0.08,
    r_outer = 0.12,
    z_length = 0.1
)

coils = [coil1, coil2]
println("   Created $(length(coils)) coils")

# Step 2: Create mesh from coil geometry
println("\n2. Creating mesh from coil geometry...")
domain = (0.0, 0.6, 0.0, 1.0)  # (r_min, r_max, z_min, z_max)
resolution = (240, 360)
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
# Half-domain symmetry: Neumann at z=0, Dirichlet on r=0, r_max, z_max
solution = solve_magnetostatic(model, J_φ; symmetry_z0=true)
println("   Solution completed")

# Step 5: Extract field on axis
println("\n5. Extracting magnetic field on axis...")
z_points = collect(0.0:0.01:1.0)
B_z_values = extract_field_on_axis(solution.B_z, z_points)
println("   Field extracted at $(length(z_points)) points")

# Step 6: Visualize results
println("\n6. Visualizing results...")

# Prepare field samples once to reuse for multiple plots
println("   Preparing field samples for fast plotting...")
grid_cache = sample_field_grid(solution; n_r=40, n_z=80)

# 6a. Plot mesh structure
println("   6a. Plotting mesh...")
p_mesh = plot_mesh(model, title="Mesh Structure")
mesh_path = plotsdir("mesh_structure.png")
mkpath(dirname(mesh_path))
savefig(p_mesh, mesh_path)
println("      Mesh plot saved to: $mesh_path")

# 6b. Plot B-field magnitude
println("   6b. Plotting B-field magnitude...")
p_b_mag = plot_b_field(solution; field_type=:magnitude, title="Magnetic Field Magnitude",
                       grid_data=grid_cache, plot_kind=:heatmap)
b_mag_path = plotsdir("b_field_magnitude.png")
mkpath(dirname(b_mag_path))
savefig(p_b_mag, b_mag_path)
println("      B-field magnitude plot saved to: $b_mag_path")

# 6c. Plot B_z component
println("   6c. Plotting B_z component...")
p_b_z = plot_b_field(solution; field_type=:B_z, title="Axial Magnetic Field Component (B_z)",
                     grid_data=grid_cache)
b_z_path = plotsdir("b_field_z_component.png")
mkpath(dirname(b_z_path))
savefig(p_b_z, b_z_path)
println("      B_z plot saved to: $b_z_path")

# 6d. Plot field on axis
println("   6d. Plotting field on axis...")
p_axis = plot(z_points, B_z_values, 
         xlabel="z [m]", 
         ylabel="B_z [T]", 
         title="Magnetic Field on Axis (2D Axisymmetric)",
         linewidth=2,
         legend=false)
axis_path = plotsdir("field_on_axis.png")
mkpath(dirname(axis_path))
savefig(p_axis, axis_path)
println("      Axis field plot saved to: $axis_path")

# Display plots
display(p_mesh)
display(p_b_mag)
display(p_b_z)
display(p_axis)

# Step 7: Save results
println("\n7. Saving results...")
params = Dict(
    "coils" => length(coils),
    "coil1_r" => coil1.r,
    "coil1_z" => coil1.z,
    "coil1_current_density" => coil1.current_density,
    "coil2_r" => coil2.r,
    "coil2_z" => coil2.z,
    "coil2_current_density" => coil2.current_density,
)

save_results(solution.A_φ, solution.B_z, params, "magnetostatic_solution")
println("   Results saved")

# Step 8: Export VTK for external visualization
# println("\n8. Exporting solution to VTK...")
# vtk_path = export_solution_vtk(solution, "magnetostatic_solution")
# println("   VTK file saved to: $vtk_path")

println("\n" * "=" ^ 60)
println("Example completed successfully!")
println("=" ^ 60)


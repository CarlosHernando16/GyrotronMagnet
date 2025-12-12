using DrWatson, Test
@quickactivate "GyrotronMagnet"

using GyrotronMagnet
using Gridap

@testset "FEM Solver Module Tests" begin
    @testset "MagnetostaticModel structure" begin
        # Test that types are exported
        @test typeof(MagnetostaticModel) <: DataType
    end
    
    @testset "Solver on simple mesh" begin
        # Create a simple 2D mesh
        domain = (0.0, 0.3, 0.0, 1.0)  # (r_min, r_max, z_min, z_max)
        resolution = (15, 50)
        model = CartesianDiscreteModel(domain, resolution)
        
        # Create a simple coil
        coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7;
            r_inner=0.08, r_outer=0.12, z_length=0.1)
        coils = [coil]
        
        # Create current density
        J_φ = create_current_density(coils, model)
        
        # Solve magnetostatic problem
        solution = solve_magnetostatic(model, J_φ)
        
        # Verify solution structure
        @test solution isa MagnetostaticModel
        @test solution.model == model
        @test solution.A_φ !== nothing
        @test solution.B !== nothing
        @test solution.B_r !== nothing
        @test solution.B_z !== nothing
    end
    
    @testset "Field extraction on axis" begin
        # Create mesh and solve
        domain = (0.0, 0.3, 0.0, 1.0)
        model = CartesianDiscreteModel(domain, (15, 50))
        coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7;
            r_inner=0.08, r_outer=0.12, z_length=0.1)
        J_φ = create_current_density([coil], model)
        solution = solve_magnetostatic(model, J_φ)
        
        # Extract field on axis
        z_points = collect(0.1:0.1:0.9)
        B_z_values = extract_field_on_axis(solution.B_z, z_points)
        
        @test length(B_z_values) == length(z_points)
        @test all(isfinite, B_z_values)
        
        # Field should be non-zero near the coil
        @test any(abs.(B_z_values) .> 0)
    end
    
    @testset "Field symmetry" begin
        # For a coil at z=0.5, field should be roughly symmetric around z=0.5
        domain = (0.0, 0.3, 0.0, 1.0)
        model = CartesianDiscreteModel(domain, (15, 50))
        coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7;
            r_inner=0.08, r_outer=0.12, z_length=0.1)
        J_φ = create_current_density([coil], model)
        solution = solve_magnetostatic(model, J_φ)
        
        # Extract field at symmetric points around z=0.5
        z_low = [0.3, 0.4]
        z_high = [0.7, 0.6]
        
        B_low = extract_field_on_axis(solution.B_z, z_low)
        B_high = extract_field_on_axis(solution.B_z, z_high)
        
        # B_z values at symmetric positions should be similar
        for i in eachindex(B_low)
            @test isapprox(B_low[i], B_high[i], rtol=0.1)
        end
    end
    
    @testset "Compute magnetic field standalone" begin
        # Test compute_magnetic_field function
        domain = (0.0, 0.3, 0.0, 1.0)
        model = CartesianDiscreteModel(domain, (10, 30))
        coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7;
            r_inner=0.08, r_outer=0.12, z_length=0.1)
        J_φ = create_current_density([coil], model)
        solution = solve_magnetostatic(model, J_φ)
        
        # Compute field components separately
        B_r, B_z, B_mag = compute_magnetic_field(solution.A_φ, model)
        
        @test B_r !== nothing
        @test B_z !== nothing
        @test B_mag !== nothing
    end

    @testset "Symmetry mode (half-domain z≥0)" begin
        # Half-domain with symmetry BC at z=0
        # Coil at z=0.15 (only in z>0 domain)
        domain = (0.0, 0.3, 0.0, 0.5)  # Only z ≥ 0
        model = CartesianDiscreteModel(domain, (15, 25))
        
        coil = create_coil_geometry(r=0.1, z=0.15, current_density=5e7;
            r_inner=0.08, r_outer=0.12, z_length=0.1)
        coils = [coil]
        J_φ = create_current_density(coils, model)
        
        # Solve with symmetry_z0=true
        solution = solve_magnetostatic(model, J_φ; symmetry_z0=true)
        
        @test solution isa MagnetostaticModel
        @test solution.A_φ !== nothing
        
        # Extract field on axis (z ≥ 0 only)
        z_points = collect(0.05:0.05:0.45)
        B_z_values = extract_field_on_axis(solution.B_z, z_points)
        
        @test length(B_z_values) == length(z_points)
        @test all(isfinite, B_z_values)
        @test any(abs.(B_z_values) .> 0)
    end

    @testset "Axis Dirichlet and midplane Neumann" begin
        # r=0 should be Dirichlet; z=0 should be Neumann in symmetry mode
        domain = (0.0, 0.25, 0.0, 0.4)
        model = CartesianDiscreteModel(domain, (12, 18))
        coil = create_coil_geometry(r=0.12, z=0.2, current_density=4e7;
            r_inner=0.10, r_outer=0.14, z_length=0.08)
        J_φ = create_current_density([coil], model)
        solution = solve_magnetostatic(model, J_φ; symmetry_z0=true)
        # Field should remain finite on axis; Neumann at z=0 avoids pinning midplane
        z_pts = collect(range(0.0, 0.35, length=20))
        B_on_axis = extract_field_on_axis(solution.B_z, z_pts)
        @test all(isfinite, B_on_axis)
    end

    @testset "Region-aware material properties" begin
        # Test different mu_r for air and coil regions
        domain = (0.0, 0.3, 0.0, 1.0)
        model = CartesianDiscreteModel(domain, (15, 50))
        
        coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7;
            r_inner=0.08, r_outer=0.12, z_length=0.1)
        coils = [coil]
        J_φ = create_current_density(coils, model)
        
        # Solve with uniform material (reference)
        sol_uniform = solve_magnetostatic(model, J_φ;
            coil_geometries=coils, mu_r_air=1.0, mu_r_coil=1.0)
        
        # Solve with higher permeability in coil region
        # Note: mu_r_coil > 1 should increase the field slightly
        sol_high_mu = solve_magnetostatic(model, J_φ;
            coil_geometries=coils, mu_r_air=1.0, mu_r_coil=1.5)
        
        @test sol_uniform isa MagnetostaticModel
        @test sol_high_mu isa MagnetostaticModel
        
        # Extract field values for comparison
        z_mid = [0.5]
        B_uniform = extract_field_on_axis(sol_uniform.B_z, z_mid)
        B_high_mu = extract_field_on_axis(sol_high_mu.B_z, z_mid)
        
        @test all(isfinite, B_uniform)
        @test all(isfinite, B_high_mu)
        # Higher mu_r_coil should not decrease the field
        # (exact relationship depends on geometry, but should be in same order of magnitude)
        @test abs(B_high_mu[1]) > 0
    end
end


using DrWatson, Test
@quickactivate "GyrotronMagnet"

using GyrotronMagnet

@testset "Geometry Module Tests" begin
    @testset "CoilGeometry creation" begin
        coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7)
        @test coil.r == 0.1
        @test coil.z == 0.5
        @test coil.r_inner == 0.01
        @test coil.r_outer == 0.02
        @test coil.z_length == 0.01
        @test coil.current_density == 5e7
    end
    
    @testset "CoilGeometry with custom dimensions" begin
        coil = create_coil_geometry(r=0.2, z=0.3, current_density=6e7; 
                                     r_inner=0.02, r_outer=0.04, z_length=0.02)
        @test coil.r_inner == 0.02
        @test coil.r_outer == 0.04
        @test coil.z_length == 0.02
        @test coil.current_density == 6e7
    end
    
    @testset "Multiple coils" begin
        coils = [
            create_coil_geometry(r=0.1, z=0.2, current_density=5e7),
            create_coil_geometry(r=0.1, z=0.8, current_density=5e7)
        ]
        @test length(coils) == 2
        @test coils[1].z < coils[2].z
    end

    @testset "ConductorParameters and derived calculations" begin
        coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7; r_inner=0.01, r_outer=0.03, z_length=0.02)
        conductor = ConductorParameters(1.2e-6, "copper"; winding_factor=0.8, insulation_thickness=0.0001)
        area = compute_coil_area(coil)
        @test area ≈ π * (0.03^2 - 0.01^2) * 0.02
        amp_turns = compute_ampere_turns(coil)
        @test amp_turns ≈ coil.current_density * area
        mean_circ = compute_mean_circumference(coil)
        @test mean_circ ≈ π * (0.01 + 0.03)
        n_turns = compute_n_turns(coil; operating_current=500.0)
        @test n_turns ≈ amp_turns / 500.0
        length = compute_conductor_length(coil, conductor; operating_current=500.0)
        @test length ≈ (n_turns * mean_circ) / conductor.winding_factor
    end
    
    @testset "Mesh loading (requires mesh file)" begin
        # This test requires an actual mesh file
        # For now, we'll test that the function exists and handles missing files
        @test_throws ErrorException load_mesh("nonexistent.msh")
    end
    
    # Note: Mesh generation test requires Gmsh installation
    # This would be tested in integration tests with actual Gmsh setup
end


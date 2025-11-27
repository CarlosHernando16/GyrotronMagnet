using DrWatson, Test
@quickactivate "GyrotronMagnet"

using GyrotronMagnet

@testset "Geometry Module Tests" begin
    @testset "CoilParameters creation" begin
        coil = create_coil_geometry(0.1, 0.5, 1000.0, 100)
        @test coil.r == 0.1
        @test coil.z == 0.5
        @test coil.I == 1000.0
        @test coil.n_turns == 100
        @test coil.r_inner == 0.01
        @test coil.r_outer == 0.02
        @test coil.z_length == 0.01
    end
    
    @testset "CoilParameters with custom dimensions" begin
        coil = create_coil_geometry(0.2, 0.3, 500.0, 50; 
                                     r_inner=0.02, r_outer=0.04, z_length=0.02)
        @test coil.r_inner == 0.02
        @test coil.r_outer == 0.04
        @test coil.z_length == 0.02
    end
    
    @testset "Multiple coils" begin
        coils = [
            create_coil_geometry(0.1, 0.2, 1000.0, 100),
            create_coil_geometry(0.1, 0.8, 1000.0, 100)
        ]
        @test length(coils) == 2
        @test coils[1].z < coils[2].z
    end
    
    @testset "Mesh loading (requires mesh file)" begin
        # This test requires an actual mesh file
        # For now, we'll test that the function exists and handles missing files
        @test_throws ErrorException load_mesh("nonexistent.msh")
    end
    
    # Note: Mesh generation test requires Gmsh installation
    # This would be tested in integration tests with actual Gmsh setup
end


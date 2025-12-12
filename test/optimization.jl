using Test
using GyrotronMagnet
using Gridap

@testset "Optimization Module Tests" begin
    @testset "Requirements" begin
        req = default_requirements()
        @test req isa OptimizationRequirements
        @test validate_requirements(req) == true
        
        # Test target profile
        profile = default_target_profile()
        @test length(profile.z) == length(profile.Bz)
        @test issorted(profile.z)
        @test profile.symmetry == true
    end
    
    @testset "Objective Function" begin
        req = default_requirements()
        
        # Create mock metrics
        z_samples = collect(range(-0.3, 0.3, length=50))
        Bz_samples = fill(2.0, length(z_samples))
        metrics = SimulationMetrics(z_samples, Bz_samples;
                                   max_conductor_field=5.0,
                                   max_radial_field=2.0)
        
        # Create mock coils
        coils = [
            CoilDescriptor(0.1, 0.15, 0.05, 0.0, 5.0e7),
            CoilDescriptor(0.1, 0.15, 0.05, 0.1, 5.0e7)
        ]
        
        total, base_mse, penalties = objective_function(metrics, coils, req)
        @test total >= 0
        @test base_mse >= 0
        @test penalties isa PenaltyBreakdown
    end
    
    @testset "Parameter Bounds" begin
        bounds = CoilParameterBounds(3;
                                    r_inner=(0.05, 0.2),
                                    r_outer=(0.08, 0.35),
                                    height=(0.01, 0.08),
                                    z_center=(-0.25, 0.25),
                                    current_density=(2.0e7, 8.0e7))
        @test length(bounds.r_inner) == 3
        @test length(bounds.r_outer) == 3
        @test length(bounds.height) == 3
        @test length(bounds.z_center) == 3
        @test length(bounds.current_density) == 3
    end
    
    @testset "FEM Evaluator" begin
        # Create a simple mesh covering the z range from default profile
        # default_target_profile uses z_max=0.5, so symmetric profile goes from -0.5 to 0.5
        domain = (0.0, 0.5, -0.6, 0.6)  # (r_min, r_max, z_min, z_max) - cover symmetric profile
        model = CartesianDiscreteModel(domain, (10, 30))
        req = default_requirements()
        
        conductor = ConductorParameters(1e-6, "copper"; winding_factor=0.7, insulation_thickness=0.0)
        evaluator = build_fem_evaluator(model, req; conductor=conductor)
        @test evaluator isa FEMEvaluator
        @test evaluator.model == model
        @test length(evaluator.z_samples) > 0
        
        # Test evaluator call
        coils = [
            CoilDescriptor(0.1, 0.15, 0.05, 0.0, 4.0e7),  # current_density in A/m²
            CoilDescriptor(0.1, 0.15, 0.05, 0.1, 4.0e7)
        ]
        
        metrics = evaluator(coils)
        @test metrics isa SimulationMetrics
        @test length(metrics.z_samples) == length(metrics.Bz_samples)
    end
    
    @testset "Optimization Config" begin
        config = OptimizationConfig(algorithm=:de, population=20, iterations=50)
        @test config.algorithm == :de
        @test config.population == 20
        @test config.iterations == 50
    end

    @testset "FEM Evaluator with symmetry mode" begin
        # Half-domain (z >= 0) with symmetry at z=0
        domain = (0.0, 0.5, 0.0, 0.6)  # Only z >= 0
        model = CartesianDiscreteModel(domain, (10, 20))
        
        # Create custom requirements with only z >= 0 target profile
        z_samples = collect(range(0.0, 0.5, length=25))
        Bz_target = fill(2.0, length(z_samples))
        profile = TargetFieldProfile(z_samples, Bz_target; symmetry=false)
        
        req = OptimizationRequirements(
            profile,
            default_requirements().em,
            default_requirements().geometry,
            default_requirements().safety,
            default_requirements().penalties
        )
        
        conductor = ConductorParameters(1e-6, "copper")
        evaluator = build_fem_evaluator(model, req;
            conductor=conductor,
            symmetry_z0=true)  # Enable symmetry mode
        
        @test evaluator isa FEMEvaluator
        @test evaluator.symmetry_z0 == true
        
        # Test evaluator call with coils in z > 0
        coils = [
            CoilDescriptor(0.1, 0.15, 0.05, 0.15, 4.0e7),  # z_center > 0
            CoilDescriptor(0.1, 0.15, 0.05, 0.35, 4.0e7)
        ]
        
        metrics = evaluator(coils)
        @test metrics isa SimulationMetrics
        @test all(isfinite, metrics.Bz_samples)
    end

    @testset "FEM Evaluator with material properties" begin
        domain = (0.0, 0.5, -0.6, 0.6)
        model = CartesianDiscreteModel(domain, (10, 30))
        req = default_requirements()
        
        conductor = ConductorParameters(1e-6, "copper")
        evaluator = build_fem_evaluator(model, req;
            conductor=conductor,
            mu_r_air=1.0,
            mu_r_coil=1.2)  # Slightly higher mu_r in coil
        
        @test evaluator isa FEMEvaluator
        @test evaluator.mu_r_air == 1.0
        @test evaluator.mu_r_coil == 1.2
        
        coils = [
            CoilDescriptor(0.1, 0.15, 0.05, 0.0, 4.0e7),
            CoilDescriptor(0.1, 0.15, 0.05, 0.1, 4.0e7)
        ]
        
        metrics = evaluator(coils)
        @test metrics isa SimulationMetrics
        @test all(isfinite, metrics.Bz_samples)
    end
    
    # Integration test skipped - requires long runtime
    # Uncomment to run full optimization test:
    # @testset "Integration: Quick Optimization" begin
    #     domain = (0.0, 0.5, -0.6, 0.6)
    #     model = CartesianDiscreteModel(domain, (8, 24))
    #     req = default_requirements()
    #     bounds = CoilParameterBounds(2;
    #                                 r_inner=(0.08, 0.12),
    #                                 r_outer=(0.12, 0.18),
    #                                 height=(0.03, 0.07),
    #                                 z_center=(-0.1, 0.1),
    #                                 current_density=(2.0e7, 6.0e7))
    #     evaluator = build_fem_evaluator(model, req; turns_per_coil=10)
    #     config = OptimizationConfig(algorithm=:de, population=5, iterations=3, verbose=false)
    #     result = run_metaheuristic(req, bounds, coils -> evaluator(coils); config=config)
    #     @test result isa OptimizationResult
    #     @test result.best_value >= 0
    #     @test length(result.coils) == 2
    # end
end


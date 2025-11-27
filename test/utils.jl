using DrWatson, Test
@quickactivate "GyrotronMagnet"

using GyrotronMagnet

@testset "Utils Module Tests" begin
    @testset "RMS error computation" begin
        B_computed = [1.0, 2.0, 3.0, 4.0]
        B_target = [1.1, 2.1, 2.9, 4.1]
        
        rms = compute_rms_error(B_computed, B_target)
        @test rms ≈ sqrt((0.1^2 + 0.1^2 + 0.1^2 + 0.1^2) / 4) atol=1e-10
        
        # Perfect match
        rms_zero = compute_rms_error(B_computed, B_computed)
        @test rms_zero == 0.0
    end
    
    @testset "RMS error with mismatched lengths" begin
        B_computed = [1.0, 2.0]
        B_target = [1.0, 2.0, 3.0]
        @test_throws ErrorException compute_rms_error(B_computed, B_target)
    end
    
    @testset "Load target profile (requires CSV file)" begin
        # Test error handling for missing file
        @test_throws ErrorException load_target_profile("nonexistent.csv")
    end
    
    # Note: Field extraction tests require actual B field from solver
    # These would be integration tests with solved magnetostatic problem
end


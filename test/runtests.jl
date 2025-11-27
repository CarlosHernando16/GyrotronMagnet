using DrWatson, Test
@quickactivate "GyrotronMagnet"

# Include test files
include("geometry.jl")
include("fem_solver.jl")
include("utils.jl")

# Run test suite
println("Starting tests")
ti = time()

@testset "GyrotronMagnet tests" begin
    # All tests are in individual test files
    # This testset serves as a container
    @test true
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")

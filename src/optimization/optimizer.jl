module OptimizationEngine

using Metaheuristics

# Import from sibling modules
using ..OptimizationSpecs: OptimizationRequirements
using ..OptimizationObjective: CoilDescriptor, SimulationMetrics, PenaltyBreakdown,
    objective_function

# Import items from parent module (GyrotronMagnet)
# These are defined in geometry.jl, fem_solver.jl, utils.jl which are included
# directly into the GyrotronMagnet module
import ..CoilGeometry
import ..ConductorParameters
import ..create_current_density
import ..solve_magnetostatic
import ..extract_field_on_axis
import ..sample_field_grid

export CoilParameterBounds, OptimizationConfig, OptimizationResult,
       FEMEvaluator, build_fem_evaluator, run_metaheuristic

"""
    CoilParameterBounds

Per-coil parameter ranges defining the optimizer search space.
All tuples are `(min, max)` in meters (or A/m² for `current_density`).
"""
struct CoilParameterBounds
    r_inner::Vector{Tuple{Float64, Float64}}
    r_outer::Vector{Tuple{Float64, Float64}}
    height::Vector{Tuple{Float64, Float64}}
    z_center::Vector{Tuple{Float64, Float64}}
    current_density::Vector{Tuple{Float64, Float64}}  # Current density bounds [A/m²]
    min_radial_gap::Float64
end

function CoilParameterBounds(n_coils::Integer;
                             r_inner::Union{Tuple,Vector{Tuple{Float64,Float64}}}=(0.05, 0.2),
                             r_outer::Union{Tuple,Vector{Tuple{Float64,Float64}}}=(0.08, 0.35),
                             height::Union{Tuple,Vector{Tuple{Float64,Float64}}}=(0.01, 0.08),
                             z_center::Union{Tuple,Vector{Tuple{Float64,Float64}}}=(-0.25, 0.25),
                             current_density::Union{Tuple,Vector{Tuple{Float64,Float64}}}=(2.0e7, 8.0e7),
                             min_radial_gap::Float64=2e-3)
    n = Int(n_coils)
    n > 0 || error("Number of coils must be positive.")
    return CoilParameterBounds(
        _spread_bounds(r_inner, n),
        _spread_bounds(r_outer, n),
        _spread_bounds(height, n),
        _spread_bounds(z_center, n),
        _spread_bounds(current_density, n),
        min_radial_gap
    )
end

"""
    OptimizationConfig(; algorithm=:de, population=0, iterations=500,
                         seed=0x12345678, verbose=false, F=0.7, CR=0.5)

High-level configuration for the metaheuristic search.
"""
Base.@kwdef struct OptimizationConfig
    algorithm::Symbol = :de
    population::Int = 0
    iterations::Int = 500
    seed::UInt = 0x12345678
    verbose::Bool = false
    F::Float64 = 0.7
    CR::Float64 = 0.6
end

"""
    OptimizationResult

Holds the best solution returned by the optimizer along with diagnostics.
"""
struct OptimizationResult
    best_value::Float64
    base_mse::Float64
    penalties::PenaltyBreakdown
    coils::Vector{CoilDescriptor}
    metrics::SimulationMetrics
    iterations::Int
    f_calls::Int
    stop_code::Int
end

"""
    FEMEvaluator(model, config)

Callable wrapper that converts a vector of `CoilDescriptor`s into a FEM solve
and returns `SimulationMetrics`.

Supports half-domain symmetry (`symmetry_z0=true`) and region-aware material properties
(`mu_r_air`, `mu_r_coil`).
"""
struct FEMEvaluator
    model
    z_samples::Vector{Float64}
    μ₀::Float64
    order::Int
    bc_type::String
    boundary_tag::String
    conductor::ConductorParameters
    sampling_grid::Tuple{Int, Int}
    r_limits::Union{Nothing,Tuple{Float64,Float64}}
    z_limits::Union{Nothing,Tuple{Float64,Float64}}
    symmetry_z0::Bool
    mu_r_air::Float64
    mu_r_coil::Float64
end

"""
    build_fem_evaluator(model; z_samples, μ₀, order, bc_type, boundary_tag,
                        conductor, sampling_grid, r_limits, z_limits,
                        symmetry_z0, mu_r_air, mu_r_coil)

Convenience constructor for `FEMEvaluator`.

# New Parameters
- `symmetry_z0::Bool`: If true, solve half-domain (z≥0) with symmetry BC at z=0 (default false)
- `mu_r_air::Float64`: Relative permeability of air region (default 1.0)
- `mu_r_coil::Float64`: Relative permeability of coil region (default 1.0)
"""
function build_fem_evaluator(model;
                             z_samples::AbstractVector{<:Real},
                             μ₀::Float64 = 4π * 1e-7,
                             order::Int = 1,
                             bc_type::String = "dirichlet",
                             boundary_tag::String = "outer_boundary",
                             conductor::ConductorParameters = ConductorParameters(1e-6, "copper"; winding_factor=0.7, insulation_thickness=0.0),
                             sampling_grid::Tuple{Int,Int} = (40, 80),
                             r_limits=nothing,
                             z_limits=nothing,
                             symmetry_z0::Bool=false,
                             mu_r_air::Float64=1.0,
                             mu_r_coil::Float64=1.0)
    z_vec = collect(float.(z_samples))
    return FEMEvaluator(model, z_vec, μ₀, order, bc_type, boundary_tag,
                        conductor, sampling_grid, r_limits, z_limits,
                        symmetry_z0, mu_r_air, mu_r_coil)
end

"""
    build_fem_evaluator(model, req::OptimizationRequirements; ...)

Convenience constructor that extracts z_samples from OptimizationRequirements.
"""
function build_fem_evaluator(model, req::OptimizationRequirements;
                             conductor::ConductorParameters = ConductorParameters(1e-6, "copper"; winding_factor=0.7, insulation_thickness=0.0),
                             μ₀::Float64 = 4π * 1e-7,
                             order::Int = 1,
                             bc_type::String = "dirichlet",
                             boundary_tag::String = "outer_boundary",
                             sampling_grid::Tuple{Int,Int} = (40, 80),
                             r_limits=nothing,
                             z_limits=nothing,
                             symmetry_z0::Bool=false,
                             mu_r_air::Float64=1.0,
                             mu_r_coil::Float64=1.0)
    return build_fem_evaluator(model; z_samples=req.target_profile.z,
                               μ₀=μ₀, order=order, bc_type=bc_type,
                               boundary_tag=boundary_tag,
                               conductor=conductor,
                               sampling_grid=sampling_grid,
                               r_limits=r_limits, z_limits=z_limits,
                               symmetry_z0=symmetry_z0,
                               mu_r_air=mu_r_air, mu_r_coil=mu_r_coil)
end

function (eval::FEMEvaluator)(coils::Vector{CoilDescriptor})
    isempty(coils) && error("At least one coil descriptor is required.")
    coil_geoms = map(_descriptor_to_geometry, coils)
    Jφ = create_current_density(coil_geoms, eval.model)
    solution = solve_magnetostatic(eval.model, Jφ;
                                   μ₀=eval.μ₀,
                                   order=eval.order,
                                   bc_type=eval.bc_type,
                                   boundary_tag=eval.boundary_tag,
                                   symmetry_z0=eval.symmetry_z0,
                                   coil_geometries=coil_geoms,
                                   mu_r_air=eval.mu_r_air,
                                   mu_r_coil=eval.mu_r_coil)
    Bz_axis = extract_field_on_axis(solution.B_z, eval.z_samples)
    grid = sample_field_grid(solution;
                             n_r=eval.sampling_grid[1],
                             n_z=eval.sampling_grid[2],
                             r_limits=eval.r_limits,
                             z_limits=eval.z_limits)
    max_mag = maximum(grid.data[:magnitude])
    max_radial = maximum(abs, grid.data[:B_r])
    return SimulationMetrics(eval.z_samples, Bz_axis;
                             max_conductor_field=max_mag,
                             max_radial_field=max_radial)
end

# -----------------------------------------------------------------------------#
# Metaheuristic driver

"""
    run_metaheuristic(req, bounds, evaluator; config=OptimizationConfig(), logger=nothing)

Optimize coil descriptors to satisfy `req` using a metaheuristic search.
"""
function run_metaheuristic(req::OptimizationRequirements,
                           bounds::CoilParameterBounds,
                           evaluator;
                           config::OptimizationConfig=OptimizationConfig(),
                           logger::Function = (status) -> nothing)
    lower, upper = _flatten_bounds(bounds)
    space = _bounds_matrix(lower, upper)
    algorithm = _build_algorithm(config, length(lower))

    fitness = let bounds=bounds, req=req, evaluator=evaluator
        function (x)
            coils = _decode_candidate(x, bounds)
            local metrics
            try
                metrics = evaluator(coils)
            catch err
                @warn "Evaluator failed, returning Inf cost" exception=(err, catch_backtrace())
                return Inf
            end
            total, _, _ = objective_function(metrics, coils, req)
            return total
        end
    end

    status = optimize(fitness, space, algorithm; logger=logger)
    best_x = status.best_sol.x
    best_coils = _decode_candidate(best_x, bounds)
    metrics = evaluator(best_coils)
    best_value, base_mse, penalties = objective_function(metrics, best_coils, req)
    # Convert termination status to integer code
    stop_code_int = if status.termination_status_code isa Integer
        Int(status.termination_status_code)
    elseif hasfield(typeof(status.termination_status_code), :code)
        Int(status.termination_status_code.code)
    else
        # Default to iteration limit if can't determine
        0
    end

    return OptimizationResult(best_value, base_mse, penalties, best_coils, metrics,
                              status.iteration, status.f_calls, stop_code_int)
end

# -----------------------------------------------------------------------------#
# Internal helpers

function _spread_bounds(values, n::Int)
    if values isa Vector
        length(values) == n ||
            error("Vector bounds must match number of coils ($n).")
        return Vector{Tuple{Float64,Float64}}(values)
    elseif values isa Tuple
        lo, hi = float(values[1]), float(values[2])
        return fill((lo, hi), n)
    else
        error("Bounds must be a tuple or vector of tuples.")
    end
end

function _flatten_bounds(bounds::CoilParameterBounds)
    lower = Float64[]
    upper = Float64[]
    for i in eachindex(bounds.r_inner)
        ri = bounds.r_inner[i]
        ro = bounds.r_outer[i]
        h = bounds.height[i]
        zc = bounds.z_center[i]
        jd = bounds.current_density[i]
        append!(lower, (ri[1], ro[1], h[1], zc[1], jd[1]))
        append!(upper, (ri[2], ro[2], h[2], zc[2], jd[2]))
    end
    return lower, upper
end

function _bounds_matrix(lower::Vector{Float64}, upper::Vector{Float64})
    length(lower) == length(upper) || error("Lower/upper bounds dimension mismatch.")
    d = length(lower)
    mat = zeros(2, d)
    @inbounds mat[1, :] = lower
    @inbounds mat[2, :] = upper
    return mat
end

function _decode_candidate(x::AbstractVector{<:Real}, bounds::CoilParameterBounds)
    n = length(bounds.r_inner)
    expected = 5n
    length(x) == expected || error("Candidate dimension $(length(x)) does not match bounds ($expected).")
    coils = Vector{CoilDescriptor}(undef, n)
    idx = 1
    for i in 1:n
        r_in = clamp(float(x[idx]), bounds.r_inner[i][1], bounds.r_inner[i][2]); idx += 1
        r_out = clamp(float(x[idx]), bounds.r_outer[i][1], bounds.r_outer[i][2]); idx += 1
        r_out = max(r_out, r_in + bounds.min_radial_gap)
        h = clamp(float(x[idx]), bounds.height[i][1], bounds.height[i][2]); idx += 1
        zc = clamp(float(x[idx]), bounds.z_center[i][1], bounds.z_center[i][2]); idx += 1
        jd = clamp(float(x[idx]), bounds.current_density[i][1], bounds.current_density[i][2]); idx += 1
        coils[i] = CoilDescriptor(r_in, r_out, h, zc, jd)
    end
    return coils
end

function _descriptor_to_geometry(coil::CoilDescriptor)
    radial_center = (coil.r_inner + coil.r_outer) / 2
    return CoilGeometry(radial_center, coil.z_center,
                        coil.r_inner, coil.r_outer, coil.height,
                        coil.current_density)
end

function _build_algorithm(config::OptimizationConfig, dim::Int)
    opts = Metaheuristics.Options(iterations=config.iterations,
                                  verbose=config.verbose,
                                  seed=config.seed)
    if config.algorithm == :de
        return Metaheuristics.DE(N=config.population,
                                 F=config.F,
                                 CR=config.CR,
                                 options=opts)
    else
        error("Unsupported algorithm $(config.algorithm).")
    end
end

end # module


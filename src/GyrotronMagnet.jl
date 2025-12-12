"""
    GyrotronMagnet

Electromagnetic design and optimization of gyrotron magnets using Gridap.jl.

This package provides tools for:
- Coil geometry definition and mesh generation
- Finite element solution of magnetostatic problems
- Field extraction and analysis
- Optimization of magnet configurations (Phase 2)
"""
module GyrotronMagnet

using DrWatson
using Gridap
using GridapGmsh
using LinearAlgebra

# Submodules
include("geometry.jl")
include("fem_solver.jl")
include("utils.jl")
include("optimization/requirements.jl")
include("optimization/objective.jl")
include("optimization/optimizer.jl")
include("optimization/cli_utils.jl")

using .OptimizationSpecs
using .OptimizationObjective
using .OptimizationEngine
using .CLIUtils

# Re-export key functions and types
export
    # Geometry module
    CoilGeometry,
    ConductorParameters,
    create_coil_geometry,
    generate_mesh,
    load_mesh,
    create_current_density,
    compute_coil_area,
    compute_mean_circumference,
    compute_ampere_turns,
    compute_total_current,
    compute_n_turns,
    compute_conductor_length,
    
    # FEM Solver module
    MagnetostaticModel,
    solve_magnetostatic,
    compute_magnetic_field,
    apply_boundary_conditions,
    
    # Utils module
    extract_field_on_axis,
    compute_field_profile,
    save_results,
    load_target_profile,
    compute_rms_error,
    FieldGrid,
    sample_field_grid,
    plot_mesh,
    plot_b_field,
    export_solution_vtk,
    # Optimization Requirements
    TargetFieldProfile,
    ElectromagneticRequirements,
    GeometryConstraints,
    SafetyConstraints,
    PenaltyWeights,
    OptimizationRequirements,
    default_target_profile,
    default_requirements,
    validate_requirements,
    # Optimization Objective
    CoilDescriptor,
    SimulationMetrics,
    PenaltyBreakdown,
    objective_function,
    CoilParameterBounds,
    OptimizationConfig,
    OptimizationResult,
    FEMEvaluator,
    build_fem_evaluator,
    run_metaheuristic,

    # CLI Utilities
    create_default_config,
    merge_with_defaults,
    config_to_requirements,
    bounds_from_cfg,
    create_mesh_from_cfg,
    build_evaluator_from_cfg,
    build_opt_config,
    run_and_save

end # module


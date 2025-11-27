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

# Re-export key functions and types
export
    # Geometry module
    CoilParameters,
    create_coil_geometry,
    generate_mesh,
    load_mesh,
    create_current_density,
    
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
    compute_rms_error

end # module


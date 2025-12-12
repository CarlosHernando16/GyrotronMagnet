# GyrotronMagnet

Electromagnetic design and optimization of gyrotron magnets using Gridap.jl.

## Overview

GyrotronMagnet is a Julia package for the finite element analysis and optimization of magnetostatic fields in gyrotron magnet systems. It uses [Gridap.jl](https://gridap.github.io/Gridap.jl/) for finite element computations and [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/) for reproducible scientific workflows.

## Features

- **Geometry Module**: Define coil geometries and generate/load meshes
- **FEM Solver**: Solve magnetostatic problems using the vector potential formulation
  - Axis (r=0) Dirichlet for regularity
  - **Symmetry mode**: Half-domain solving (z ≥ 0) with Neumann BC at z=0
  - **Material regions**: Spatially varying permeability for air/coil regions
- **Field Analysis**: Extract and analyze magnetic fields
- **Visualization & Export**: Cache field grids, plot meshes/B-fields, export to VTK
- **Optimization**: Optimize coil configurations (Phase 2)

## Installation

```julia
using Pkg
Pkg.add("GyrotronMagnet")
```

## Quick Start

```julia
using GyrotronMagnet

# Define coil parameters (keyword-only for clarity)
coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7)

# Load mesh
model, labels = load_mesh("mesh.msh")

# Create current density (scalar J_φ for axisymmetric)
J_φ = create_current_density([coil], model)

# Solve magnetostatic problem
solution = solve_magnetostatic(model, J_φ)

# Extract field on axis
z_points = 0.0:0.01:1.0
B_z = extract_field_on_axis(solution.B_z, collect(z_points))

# Cache samples for fast plotting
grid_cache = sample_field_grid(solution; n_r=80, n_z=160)
plot_b_field(solution; grid_data=grid_cache, plot_kind=:heatmap)

# Export to VTK for external visualization
vtk_path = export_solution_vtk(solution, "example_solution")
```

## Documentation

See the [API Reference](@ref) for detailed function documentation.

## Phase Status

**Phase 1** ✅ includes:
- Core FEM solver for magnetostatic problems
- Geometry handling and mesh support
- Field extraction utilities
- Basic documentation

**Phase 2** ✅ includes:
- Optimization requirements and constraints specification
- Objective function with penalty terms
- Metaheuristic optimization engine (Metaheuristics.jl)
- Workflow script for automated optimization
- Comprehensive testing

## Phase 2: Optimization Example

```julia
using GyrotronMagnet
using Gridap

# Define optimization requirements
req = default_requirements()

# Create mesh
domain = (0.0, 0.5, -0.6, 0.6)
model = CartesianDiscreteModel(domain, (20, 40))

# Define parameter bounds
bounds = CoilParameterBounds(4;
    r_inner=(0.08, 0.15),
    r_outer=(0.12, 0.25),
    height=(0.02, 0.08),
    z_center=(-0.2, 0.2),
    current_density=(2.0e7, 8.0e7)  # Current density [A/m²]
)

# Build evaluator and run optimization
conductor = ConductorParameters(1e-6, "copper"; winding_factor=0.7, insulation_thickness=0.0)
evaluator = build_fem_evaluator(model, req; conductor=conductor)
config = OptimizationConfig(algorithm=:de, population=30, iterations=100)
result = run_metaheuristic(req, bounds, coils -> evaluator(coils); config=config)
```

## Symmetry Mode and Material Regions

### Half-Domain Symmetry

For symmetric problems, solve only the z ≥ 0 half-domain:

```julia
# Via solver directly
solution = solve_magnetostatic(model, J_φ; symmetry_z0=true)

# Or via config
cfg = create_default_config()
cfg["mesh"]["symmetry_z0"] = true
```

### Material Properties

Specify different permeabilities for air and coil regions:

```julia
solution = solve_magnetostatic(model, J_φ;
    coil_geometries=coils,
    mu_r_air=1.0,
    mu_r_coil=1.5
)

# Or via config
cfg["materials"]["mu_r_air"] = 1.0
cfg["materials"]["mu_r_coil"] = 1.2
```

## Design note: physics vs. engineering parameters

- `CoilGeometry` holds what the FEM needs: geometry plus azimuthal current density (A/m²).
- `ConductorParameters` holds wire/packing properties; derived helpers compute turns, conductor length, and ampere-turns as needed.

## Citation

If you use GyrotronMagnet in your research, please cite:

```bibtex
@software{gyrotronmagnet2024,
  author = {Hernando, Carlos},
  title = {GyrotronMagnet: Electromagnetic Design of Gyrotron Magnets},
  year = {2024},
  url = {https://github.com/CarlosHernando16/GyrotronMagnet}
}
```


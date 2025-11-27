# GyrotronMagnet

Electromagnetic design and optimization of gyrotron magnets using Gridap.jl.

## Overview

GyrotronMagnet is a Julia package for the finite element analysis and optimization of magnetostatic fields in gyrotron magnet systems. It uses [Gridap.jl](https://gridap.github.io/Gridap.jl/) for finite element computations and [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/) for reproducible scientific workflows.

## Features

- **Geometry Module**: Define coil geometries and generate/load meshes
- **FEM Solver**: Solve magnetostatic problems using the vector potential formulation
- **Field Analysis**: Extract and analyze magnetic fields
- **Optimization**: Optimize coil configurations (Phase 2)

## Installation

```julia
using Pkg
Pkg.add("GyrotronMagnet")
```

## Quick Start

```julia
using GyrotronMagnet

# Define coil parameters
coil = create_coil_geometry(0.1, 0.5, 1000.0, 100)

# Load mesh
model, labels = load_mesh("mesh.msh")

# Create current density (scalar J_φ for axisymmetric)
J_φ = create_current_density([coil], model)

# Solve magnetostatic problem
solution = solve_magnetostatic(model, J_φ)

# Extract field on axis
z_points = 0.0:0.01:1.0
B_z = extract_field_on_axis(solution.B_z, collect(z_points))
```

## Documentation

See the [API Reference](@ref) for detailed function documentation.

## Phase 1 Status

Phase 1 includes:
- Core FEM solver for magnetostatic problems
- Geometry handling and mesh support
- Field extraction utilities
- Basic documentation

Phase 2 will add optimization capabilities.

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


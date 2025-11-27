# GyrotronMagnet

Electromagnetic design and optimization of gyrotron magnets using [Gridap.jl](https://gridap.github.io/Gridap.jl/) and [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/).

## Overview

GyrotronMagnet is a Julia package for finite element analysis and optimization of magnetostatic fields in gyrotron magnet systems. This project uses the [Julia Language](https://julialang.org/) and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) for reproducible scientific workflows.

**Author**: Carlos Hernando

## Phase 1 Status

Phase 1 includes:
- ✅ Core FEM solver for magnetostatic problems (∇×(ν∇×A) = J)
- ✅ Geometry handling and mesh support (GridapGmsh)
- ✅ Field extraction utilities
- ✅ Basic documentation (Documenter.jl)

**Phase 2** (planned) will add optimization capabilities.

## Installation

To (locally) reproduce this project:

1. Clone or download this repository
2. Open a Julia console and do:
   ```julia
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/GyrotronMagnet")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages including Gridap, GridapGmsh, and other dependencies.

## Quick Start

### Basic Usage

```julia
using DrWatson
@quickactivate "GyrotronMagnet"
using GyrotronMagnet

# Define coil geometry (r, z, I, n_turns; r_inner, r_outer, z_length)
coil = create_coil_geometry(0.1, 0.5, 1000.0, 100;
    r_inner = 0.08,   # Inner radius [m]
    r_outer = 0.12    # Outer radius [m]
)

# Load mesh (requires .msh file)
model, labels = load_mesh("data/mesh.msh")

# Create current density field (scalar J_φ for axisymmetric)
J_φ = create_current_density([coil], model)

# Solve magnetostatic problem
solution = solve_magnetostatic(model, J_φ)

# Extract field on axis
z_points = collect(0.0:0.01:1.0)
B_z = extract_field_on_axis(solution.B_z, z_points)
```

### Example Script

Run the complete workflow example:

```julia
include(scriptsdir("example_magnetostatic.jl"))
```

## Project Structure

```
GyrotronMagnet/
├── src/
│   ├── GyrotronMagnet.jl    # Main module
│   ├── geometry.jl           # Coil geometry and mesh handling
│   ├── fem_solver.jl         # Magnetostatic FEM solver
│   └── utils.jl              # Field extraction and utilities
├── scripts/
│   └── example_magnetostatic.jl  # Example workflow
├── test/
│   ├── runtests.jl
│   ├── geometry.jl
│   ├── fem_solver.jl
│   └── utils.jl
├── docs/
│   ├── make.jl               # Documenter setup
│   └── src/
│       ├── index.md
│       └── api.md
└── Project.toml
```

## Documentation

### Building Documentation

To build the documentation locally:

```julia
julia> using Pkg
julia> Pkg.activate("GyrotronMagnet")
julia> Pkg.add("Documenter")
julia> include("docs/make.jl")
```

The documentation will be generated in `docs/build/`.

### Online Documentation

See the [API Reference](docs/src/api.md) for detailed function documentation.

## Key Features

### Geometry Module
- Define coil geometries parametrically
- Load Gmsh mesh files (.msh)
- Generate current density fields from coil parameters

### FEM Solver Module
- Solve magnetostatic problems: ∇×(ν∇×A) = J
- Weak form implementation using Gridap
- Compute magnetic field: B = ∇×A
- Boundary condition handling

### Utilities Module
- Extract field along z-axis
- Compute field profiles along arbitrary paths
- Save results using DrWatson
- Load target field profiles
- Compute RMS errors

## Testing

Run the test suite:

```julia
julia> using Pkg
julia> Pkg.test("GyrotronMagnet")
```

Or from the project directory:

```julia
include("test/runtests.jl")
```

## Physics Background

See [PHYSICS.md](PHYSICS.md) for the mathematical formulation and [ARQUITECTURE.md](ARQUITECTURE.md) for the software architecture.

## Dependencies

- **Gridap.jl**: Finite element framework
- **GridapGmsh.jl**: Gmsh mesh support
- **DrWatson.jl**: Reproducible scientific workflows
- **Documenter.jl**: Documentation generation
- **Plots.jl** / **Makie.jl**: Visualization

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

[Add your license here]

## Citation

If you use GyrotronMagnet in your research, please cite:

```bibtex
@software{gyrotronmagnet2024,
  author = {Hernando, Carlos},
  title = {GyrotronMagnet: Electromagnetic Design of Gyrotron Magnets},
  year = {2024}
}
```

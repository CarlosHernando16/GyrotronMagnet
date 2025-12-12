# GyrotronMagnet

Electromagnetic design and optimization of gyrotron magnets using [Gridap.jl](https://gridap.github.io/Gridap.jl/) and [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/).

## Overview

GyrotronMagnet is a Julia package for finite element analysis and optimization of magnetostatic fields in gyrotron magnet systems. This project uses the [Julia Language](https://julialang.org/) and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) for reproducible scientific workflows.

**Author**: Carlos Hernando

## Phase Status

**Phase 1** ✅ includes:
- Core FEM solver for magnetostatic problems (∇×(ν∇×A) = J)
- Geometry handling and mesh support (GridapGmsh)
- Field extraction utilities
- Basic documentation (Documenter.jl)

**Phase 2** ✅ includes:
- Optimization requirements and constraints specification
- Objective function with penalty terms
- Metaheuristic optimization engine (Metaheuristics.jl)
- Workflow script for automated optimization
- Comprehensive testing

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

# Define coil geometry (keyword-only for clarity)
coil = create_coil_geometry(r=0.1, z=0.5, current_density=5e7;
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

### Example Scripts

Run the complete workflow example:

```julia
include(scriptsdir("example_magnetostatic.jl"))
```

### Phase 2: Optimization

Optimize coil geometries to meet target field profiles:

```julia
using GyrotronMagnet
using Gridap

# Define optimization requirements
req = default_requirements()  # Or customize with OptimizationRequirements

# Create mesh
domain = (0.0, 0.5, -0.6, 0.6)  # (r_min, r_max, z_min, z_max)
model = CartesianDiscreteModel(domain, (20, 40))

# Define parameter bounds for coils
bounds = CoilParameterBounds(4;  # 4 coils
    r_inner=(0.08, 0.15),
    r_outer=(0.12, 0.25),
    height=(0.02, 0.08),
    z_center=(-0.2, 0.2),
    current_density=(2.0e7, 8.0e7)  # Current density [A/m²]
)

# Build FEM evaluator with conductor parameters
conductor = ConductorParameters(1e-6, "copper")
evaluator = build_fem_evaluator(model, req; conductor=conductor)

# Configure optimization
config = OptimizationConfig(algorithm=:de, population=30, iterations=100)

# Run optimization
result = run_metaheuristic(req, bounds, evaluator; config=config)

println("Best objective: $(result.best_value)")
println("Best coils: $(length(result.coils))")
```

Or use the command-line workflow script:

```bash
julia scripts/optimize_geometry.jl --quick-test
```

See `scripts/optimize_geometry.jl` for full CLI options and JSON configuration support.

### Symmetry Mode (Half-Domain)

For symmetric problems, solve only the z ≥ 0 half-domain with a symmetry boundary condition at z=0:

```julia
# Create half-domain mesh (z >= 0 only)
domain = (0.0, 0.5, 0.0, 0.6)  # r: 0-0.5m, z: 0-0.6m
model = CartesianDiscreteModel(domain, (20, 30))

# Build evaluator with symmetry enabled
evaluator = build_fem_evaluator(model, req;
    conductor=conductor,
    symmetry_z0=true  # Enable Neumann BC at z=0
)
```

Or via JSON config:
```json
{
  "mesh": {
    "domain": [[0.0, 0.5], [0.0, 0.6]],
    "resolution": 0.02,
    "symmetry_z0": true
  }
}
```

### Material Regions

Specify different relative permeabilities for air and coil regions:

```julia
# Solve with higher permeability in coil region
solution = solve_magnetostatic(model, J_φ;
    coil_geometries=coils,
    mu_r_air=1.0,
    mu_r_coil=1.5
)

# Or via evaluator
evaluator = build_fem_evaluator(model, req;
    mu_r_air=1.0,
    mu_r_coil=1.2
)
```

Via JSON config:
```json
{
  "materials": {
    "mu_r_air": 1.0,
    "mu_r_coil": 1.2
  }
}
```

## Project Structure

```
GyrotronMagnet/
├── src/
│   ├── GyrotronMagnet.jl    # Main module
│   ├── geometry.jl           # Coil geometry and mesh handling
│   ├── fem_solver.jl         # Magnetostatic FEM solver
│   ├── utils.jl              # Field extraction and utilities
│   └── optimization/
│       ├── requirements.jl   # Optimization requirements & constraints
│       ├── objective.jl      # Objective function & penalties
│       └── optimizer.jl      # Metaheuristic optimization engine
├── scripts/
│   ├── example_magnetostatic.jl  # Phase 1 example workflow
│   └── optimize_geometry.jl     # Phase 2 optimization workflow
├── test/
│   ├── runtests.jl
│   ├── geometry.jl
│   ├── fem_solver.jl
│   ├── utils.jl
│   └── optimization.jl      # Phase 2 optimization tests
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
- Axis (r=0) Dirichlet for regularity; optional z=0 symmetry uses Neumann at midplane
- **Symmetry mode**: Half-domain solving (z ≥ 0) with Neumann BC at z=0
- **Material regions**: Spatially varying permeability (μᵣ) for air and coil regions

### Utilities Module
- Extract field along z-axis
- Compute field profiles along arbitrary paths
- Sample and cache field grids for fast plotting (`FieldGrid`, `sample_field_grid`)
- Visualize meshes and fields (`plot_mesh`, `plot_b_field`)
- Save results using DrWatson or export to VTK (`save_results`, `export_solution_vtk`)
- Load target field profiles
- Compute RMS errors

### Optimization Module (Phase 2)
- **Requirements**: Define electromagnetic targets, geometric constraints, and safety limits
- **Objective Function**: Minimize field mismatch with penalty terms for constraint violations
- **Optimization Engine**: Metaheuristic search (Differential Evolution) with configurable parameters
- **Workflow Script**: Command-line tool with JSON configuration support

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

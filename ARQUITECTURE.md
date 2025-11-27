# Project Architecture

## Core Components

### 1. Create geometry, mesh or read mesh Module (`src/geometry.jl`)
- Define coil geometries parametrically
- Create mesh based geometry using gmsh and gridapgmsh
- Load Gmsh meshes (.msh files)

### 2. FEM Solver Module (`src/fem_solver.jl`)
- Set up Gridap model and FE spaces
- Assemble magnetostatic problem: ∇×(ν∇×A) = J
- Apply boundary conditions
- Solve linear system
- Compute derived quantities (B field, field gradients)

### 3. Optimization Module (`src/optimizer.jl`)
- Define parameter vector and bounds
- Implement objective function (field profile matching)
- Compute gradients via ForwardDiff
- Interface with NLopt/Optim algorithms
- Handle constraints

### 4. Utilities (`src/utils.jl`)
- Field extraction and interpolation
- Metric calculations (RMS error, max deviation)
- Data saving with DrWatson

## Workflow

1. Load target field profile (CSV or analytical function)
2. Initialize coil parameters (positions, currents)
3. For each optimization iteration:
   - Update coil geometry
   - Generate current density J
   - Solve FEM problem for A
   - Compute B = ∇×A
   - Evaluate objective: ||B - B_target||²
4. Save optimal configuration and visualize results

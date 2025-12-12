
You are upgrading the [GyrotronMagnet](https://github.com/CarlosHernando16/GyrotronMagnet) Julia project, which uses finite element analysis to design HTS magnets for gyrotrons, by adding automated geometry optimization (Phase 2). All code should follow scientific best practices, using Gridap.jl (FEM), DrWatson.jl (reproducibility), and Julia ecosystem standards.

### 1. **Requirements and Constraints**

#### Electromagnetic Field Requirements

- The central field at the cavity region must be follow an input function


#### Geometric and Spatial Constraints

- **Bore diameter:** 
- **Total magnet height:** 
- **Outer diameter:** Must not exceed spatial envelope set by gyrotron integration.
- **No ferromagnetic materials** or magnetic sources close to the magnet coil pack.


#### Additional Constraints

- (Optional) Allow user to specify maxima for: wire usage, coil count, or current density.


### 2. **Objective Function Definition**

- Implement an objective function for geometry optimization:
    - The default objective is to minimize the squared difference between simulated and target field at prescribed axial locations:

$$
\text{Objective} = \sum_i \left[B_z(z_i; \Theta) - B_{z,\text{target}}(z_i)\right]^2
$$

where $\Theta$ = coil geometry parameters.
    - Optionally, add penalty terms for exceeding geometric or field constraints (e.g., via large multipliers).
    - Allow for alternative objectives (minimize total conductor volume, or maximize field homogeneity in given regions).


### 3. **Optimization Algorithm**

- Use a global optimization algorithm robust to local minima and mixed variables (continuous and discrete). Recommend:
    - DifferentialEvolution.jl (see thesis for rationale and parameter selection).
    - Optionally, comparison with BlackBoxOptim.jl or Optim.jl.


### 4. **Workflow Integration**

- Implement main code in src/
- Implement script in scripts/
- Allow user to specify:
    - Target field profile (as a vector of axial positions and values)
    - Geometric/physical constraints (bore, height, outer diameter, max field on conductor)
    - Material properties (critical current, etc.)
- Output results in a format compatible with Gridap.jl postprocessing, saving project state via DrWatson.jl.
- Document the workflow and provide an example with pseudo-data.


### 5. **Testing and Documentation**

- Write unit tests for the objective function and the constraint handling.
- Include at least one full optimization example with mock field/profile and geometric limits.
- Document all new functions and expose API in the module.


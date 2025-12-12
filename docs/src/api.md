# API Reference

# Module

```@docs
GyrotronMagnet
```

## Geometry Module

```@docs
GyrotronMagnet.CoilGeometry
GyrotronMagnet.ConductorParameters
GyrotronMagnet.create_coil_geometry
GyrotronMagnet.generate_mesh
GyrotronMagnet.load_mesh
GyrotronMagnet.create_current_density
GyrotronMagnet.compute_coil_area
GyrotronMagnet.compute_mean_circumference
GyrotronMagnet.compute_ampere_turns
GyrotronMagnet.compute_total_current
GyrotronMagnet.compute_n_turns
GyrotronMagnet.compute_conductor_length
GyrotronMagnet.CLIUtils.create_default_config
GyrotronMagnet.CLIUtils.merge_with_defaults
GyrotronMagnet.CLIUtils.config_to_requirements
GyrotronMagnet.CLIUtils.bounds_from_cfg
GyrotronMagnet.CLIUtils.create_mesh_from_cfg
GyrotronMagnet.CLIUtils.build_evaluator_from_cfg
GyrotronMagnet.CLIUtils.build_opt_config
GyrotronMagnet.CLIUtils.run_and_save
```

## FEM Solver Module

```@docs
GyrotronMagnet.MagnetostaticModel
GyrotronMagnet.solve_magnetostatic
GyrotronMagnet.compute_magnetic_field
GyrotronMagnet.apply_boundary_conditions
```

## Utilities Module

```@docs
GyrotronMagnet.extract_field_on_axis
GyrotronMagnet.compute_field_profile
GyrotronMagnet.save_results
GyrotronMagnet.load_target_profile
GyrotronMagnet.compute_rms_error
GyrotronMagnet.FieldGrid
GyrotronMagnet.sample_field_grid
GyrotronMagnet.plot_mesh
GyrotronMagnet.plot_b_field
GyrotronMagnet.export_solution_vtk
```

## Optimization Specifications

```@docs
GyrotronMagnet.TargetFieldProfile
GyrotronMagnet.ElectromagneticRequirements
GyrotronMagnet.GeometryConstraints
GyrotronMagnet.SafetyConstraints
GyrotronMagnet.PenaltyWeights
GyrotronMagnet.OptimizationRequirements
GyrotronMagnet.default_target_profile
GyrotronMagnet.default_requirements
GyrotronMagnet.validate_requirements
```

## Optimization Objective

```@docs
GyrotronMagnet.CoilDescriptor
GyrotronMagnet.SimulationMetrics
GyrotronMagnet.PenaltyBreakdown
GyrotronMagnet.objective_function
```

## Optimization Engine

```@docs
GyrotronMagnet.CoilParameterBounds
GyrotronMagnet.OptimizationConfig
GyrotronMagnet.OptimizationResult
GyrotronMagnet.FEMEvaluator
GyrotronMagnet.build_fem_evaluator
GyrotronMagnet.run_metaheuristic
```


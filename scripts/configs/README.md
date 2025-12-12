# Configuration Examples for GyrotronMagnet Optimization

This directory contains example JSON configuration files for the optimization scripts.

## Files

- **`config_example.json`**: Comprehensive example showing all configuration options including custom target profiles, electromagnetic requirements, geometric constraints, safety limits, penalties, bounds, optimization settings, and mesh parameters.

- **`config_simple.json`**: Minimal example showing only how to specify a custom target field profile while using defaults for all other settings.

## Usage

Use these configs with the custom optimization script:

```bash
# Run with full example config
julia --project=. scripts/custom_optimize.jl --config scripts/configs/config_example.json

# Run with simple custom profile config
julia --project=. scripts/custom_optimize.jl --config scripts/configs/config_simple.json --quick-test
```

## Configuration Sections

### target_profile
Define custom z-values and B-values for the target magnetic field profile:
```json
{
  "target_profile": {
    "z": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5],
    "Bz": [2.0, 2.0, 1.8, 1.5, 1.2, 0.8],
    "symmetry": false
  }
}
```

### Other Sections
See the main documentation or script comments for details on:
- `electromagnetic`: EM requirements (center field, homogeneity, decay)
- `geometry`: Geometric constraints (bore diameter, magnet height, etc.)
- `safety`: Field and material limits
- `penalties`: Cost function weights
- `bounds`: Optimization search space
- `optimization`: Algorithm parameters
- `mesh`: Computational domain settings

## Notes
- All sections are optional; missing sections use sensible defaults
- The script merges provided config with defaults
- Units: meters (m), tesla (T), amperes per square meter (A/m²), etc.

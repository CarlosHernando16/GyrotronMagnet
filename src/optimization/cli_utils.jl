module CLIUtils

using ..OptimizationSpecs: OptimizationRequirements, GeometryConstraints,
    TargetFieldProfile, ElectromagneticRequirements, SafetyConstraints,
    PenaltyWeights, validate_requirements
using ..OptimizationObjective: CoilDescriptor, SimulationMetrics
using ..OptimizationEngine: CoilParameterBounds, OptimizationConfig, FEMEvaluator,
    run_metaheuristic, build_fem_evaluator

using JSON
using Dates
using DrWatson
using Gridap
using GyrotronMagnet

# Public API
export create_default_config, merge_with_defaults, bounds_from_cfg, create_mesh_from_cfg,
       config_to_requirements, build_evaluator_from_cfg, build_opt_config, run_and_save

"""
    create_default_config() -> Dict

Return the default configuration dictionary with all optimization parameters.
Includes mesh settings (`symmetry_z0`), materials (`mu_r_air`, `mu_r_coil`), and all other defaults.
"""
function create_default_config()
    return Dict(
        "target_profile" => Dict(
            "z" => collect(0.0:0.01:0.5),  # Default z grid [m]
            "Bz" => fill(2.0, 51),         # Default flat 2T profile [T]
            "symmetry" => false            # No symmetry by default
        ),
        "electromagnetic" => Dict(
            "center_field" => 2.0,
            "homogeneity_ppm" => 2500.0,
            "homogeneity_window" => [0.0, 0.25],
            "decay_point" => [0.195, 0.17]
        ),
        "geometry" => Dict(
            "bore_diameter" => 0.18,
            "magnet_height_min" => 0.18,
            "magnet_height_max" => 0.20,
            "outer_diameter" => 0.60,
            "alignment_tol" => 5e-4,
            "coil_count_min" => 2,
            "coil_count_max" => 12,
            "max_current_density" => 1.0e8
        ),
        "safety" => Dict(
            "max_conductor_field" => 8.0,
            "max_radial_field" => 4.0,
            "forbid_ferromagnetics" => true
        ),
        "penalties" => Dict(
            "field" => 1.0,
            "geometry" => 1.0,
            "safety" => 1.0,
            "conductor" => 1.0
        ),
        "bounds" => Dict(
            "n_coils" => 4,
            "r_inner" => [0.06, 0.15],
            "r_outer" => [0.10, 0.25],
            "height" => [0.02, 0.06],
            "z_center" => [-0.10, 0.40],
            "current_density" => [2.0e7, 6.0e7]
        ),
        "optimization" => Dict(
            "algorithm" => "de",
            "population" => 50,
            "iterations" => 50,
            "seed" => 0x12345678,
            "verbose" => true
        ),
        "mesh" => Dict(
            "domain" => [[0.0, 0.5], [0.0, 1.0]],
            "resolution" => 0.02,
            "symmetry_z0" => false  # If true, solve only z >= 0 with symmetry BC at z=0
        ),
        "materials" => Dict(
            "mu_r_air" => 1.0,   # Relative permeability of air
            "mu_r_coil" => 1.0   # Relative permeability of coil region
        )
    )
end

""" 
    merge_with_defaults(cfg::Dict)

Merge a partial user configuration with `create_default_config()`, returning a full config dict.
Use this when a user supplies a JSON file missing some keys.
"""
function merge_with_defaults(cfg::Dict)
    defaults = create_default_config()
    return merge(defaults, cfg)
end

"""
    config_to_requirements(cfg::Dict) -> OptimizationRequirements

Build an `OptimizationRequirements` instance from a config dictionary.
Expects `target_profile`, `electromagnetic`, `geometry`, `safety`, and `penalties` sections.
"""
function config_to_requirements(cfg::Dict)
    # Target profile
    profile_dict = cfg["target_profile"]
    profile = TargetFieldProfile(
        Float64.(profile_dict["z"]),
        Float64.(profile_dict["Bz"]);
        symmetry = profile_dict["symmetry"]
    )

    # Electromagnetic requirements
    em_dict = cfg["electromagnetic"]
    em = ElectromagneticRequirements(
        em_dict["center_field"],
        em_dict["homogeneity_ppm"],
        Tuple(em_dict["homogeneity_window"]),
        Tuple(em_dict["decay_point"])
    )

    # Geometry constraints
    geom_dict = cfg["geometry"]
    geom = GeometryConstraints(
        geom_dict["bore_diameter"],
        (geom_dict["magnet_height_min"], geom_dict["magnet_height_max"]),
        geom_dict["outer_diameter"],
        geom_dict["alignment_tol"],
        (geom_dict["coil_count_min"], geom_dict["coil_count_max"]),
        get(geom_dict, "max_current_density", 8.0e7)
    )

    # Safety constraints
    safety_dict = cfg["safety"]
    safety = SafetyConstraints(
        safety_dict["max_conductor_field"],
        safety_dict["max_radial_field"],
        safety_dict["forbid_ferromagnetics"]
    )

    # Penalties
    pen_dict = cfg["penalties"]
    penalties = PenaltyWeights(
        pen_dict["field"],
        pen_dict["geometry"],
        pen_dict["safety"],
        pen_dict["conductor"]
    )

    return OptimizationRequirements(profile, em, geom, safety, penalties)
end

"""
    bounds_from_cfg(bounds_cfg::Dict) -> CoilParameterBounds

Build `CoilParameterBounds` from the `bounds` section of a config dict.
Expects keys: `n_coils`, `r_inner`, `r_outer`, `height`, `z_center`, `current_density`.
Scalar entries are promoted to fixed ranges, e.g. `r_inner = 0.08` becomes `(0.08, 0.08)`.
"""
function bounds_from_cfg(bounds_cfg::Dict)
    # Allow scalar values (fixed) or length-2 vectors/tuples (ranges)
    promote_range(x) = x isa Number ? (x, x) : Tuple(x)

    # n_coils may be a single integer or a length-1 vector; normalize to Int
    n_coils_raw = bounds_cfg["n_coils"]
    n_coils = n_coils_raw isa AbstractVector ? Int(n_coils_raw[1]) : Int(n_coils_raw)

    return CoilParameterBounds(n_coils;
        r_inner = promote_range(bounds_cfg["r_inner"]),
        r_outer = promote_range(bounds_cfg["r_outer"]),
        height  = promote_range(bounds_cfg["height"]),
        z_center = promote_range(bounds_cfg["z_center"]),
        current_density = promote_range(bounds_cfg["current_density"])
    )
end

"""
    create_mesh_from_cfg(cfg::Dict, mesh_file=nothing) -> (model, symmetry_z0::Bool)

Create a `CartesianDiscreteModel` from config mesh settings or load from `mesh_file` if provided.
Config expects `mesh.domain` as nested tuples, `mesh.resolution` as spacing, and optionally
`mesh.symmetry_z0` (default false) to enable half-domain solving with symmetry at z=0.

Returns a tuple `(model, symmetry_z0)` where `symmetry_z0` indicates whether symmetry mode is active.
"""
function create_mesh_from_cfg(cfg::Dict, mesh_file=nothing)
    mesh_cfg = cfg["mesh"]
    symmetry_z0 = get(mesh_cfg, "symmetry_z0", false)

    if mesh_file !== nothing && isfile(mesh_file)
        @info "Loading mesh from file: $mesh_file"
        model, _ = load_mesh(mesh_file)
        return (model, symmetry_z0)
    end

    domain = mesh_cfg["domain"]
    resolution = mesh_cfg["resolution"]

    r_range = Tuple(domain[1])
    z_range = Tuple(domain[2])

    # For symmetry_z0, enforce z_min = 0 and solve only z >= 0
    if symmetry_z0
        z_min = 0.0
        z_max = z_range[2]
        @info "Symmetry mode: solving half-domain z ∈ [0, $z_max] with symmetry BC at z=0"
    else
        z_min = z_range[1]
        z_max = z_range[2]
    end

    n_r = max(10, Int(ceil((r_range[2] - r_range[1]) / resolution)))
    n_z = max(10, Int(ceil((z_max - z_min) / resolution)))

    @info "Creating structured mesh: r ∈ [$(r_range[1]), $(r_range[2])], z ∈ [$z_min, $z_max], resolution=$resolution"
    domain_flat = (r_range[1], r_range[2], z_min, z_max)
    model = CartesianDiscreteModel(domain_flat, (n_r, n_z))

    return (model, symmetry_z0)
end

"""
    build_evaluator_from_cfg(model, req, cfg; symmetry_z0=false) -> FEMEvaluator

Construct an `FEMEvaluator` using conductor parameters and material properties from config.
Pass `symmetry_z0=true` to enable half-domain solving with symmetry BC at z=0.

Config sections used:
- `conductor`: wire properties (optional, defaults to copper)
- `materials`: `mu_r_air`, `mu_r_coil` (optional, default 1.0)
"""
function build_evaluator_from_cfg(model, req::OptimizationRequirements, cfg::Dict; symmetry_z0::Bool=false)
    # Conductor parameters
    conductor_cfg = get(cfg, "conductor", Dict(
        "wire_cross_section" => 1e-6,
        "material" => "copper",
        "winding_factor" => 0.7,
        "insulation_thickness" => 0.0
    ))
    conductor = ConductorParameters(
        conductor_cfg["wire_cross_section"],
        String(conductor_cfg["material"]);
        winding_factor = conductor_cfg["winding_factor"],
        insulation_thickness = conductor_cfg["insulation_thickness"]
    )

    # Material properties
    materials_cfg = get(cfg, "materials", Dict("mu_r_air" => 1.0, "mu_r_coil" => 1.0))
    mu_r_air = get(materials_cfg, "mu_r_air", 1.0)
    mu_r_coil = get(materials_cfg, "mu_r_coil", 1.0)

    return build_fem_evaluator(model, req;
        conductor=conductor,
        symmetry_z0=symmetry_z0,
        mu_r_air=mu_r_air,
        mu_r_coil=mu_r_coil
    )
end

"""
    build_opt_config(cfg::Dict) -> OptimizationConfig

Construct `OptimizationConfig` from the `optimization` section of the config dict.
"""
function build_opt_config(cfg::Dict)
    opt_cfg = cfg["optimization"]
    return OptimizationConfig(;
        algorithm = Symbol(opt_cfg["algorithm"]),
        population = opt_cfg["population"],
        iterations = opt_cfg["iterations"],
        seed = opt_cfg["seed"],
        verbose = opt_cfg["verbose"]
    )
end

"""
    run_and_save(cfg; mesh_file=nothing, output_dir=nothing, quick=false) -> NamedTuple

End-to-end pipeline: merge config, build requirements/bounds/mesh/evaluator, run optimizer, and save artifacts (config, best coils, field profile).
Returns `(result=result, out_dir=out_dir)` where `result` is `OptimizationResult`.
"""
function run_and_save(cfg::Dict; mesh_file=nothing, output_dir=nothing, quick=false)
    # Quick overrides
    if quick
        cfg["optimization"]["iterations"] = min(cfg["optimization"]["iterations"], 10)
        cfg["optimization"]["population"] = min(cfg["optimization"]["population"], 10)
        @info "Quick mode: iterations=$(cfg["optimization"]["iterations"]), population=$(cfg["optimization"]["population"])"
    end

    # Prepare output
    out_dir = output_dir === nothing ? datadir("optimization", "runs", "run_$(Dates.format(now(), "yyyy-mm-dd_HH-MM-SS"))") : output_dir
    mkpath(out_dir)

    # Save input config
    open(joinpath(out_dir, "config.json"), "w") do f
        JSON.print(f, cfg, 4)
    end

    # Requirements and bounds
    req = config_to_requirements(cfg)
    validate_requirements(req)
    bounds = bounds_from_cfg(cfg["bounds"])

    # Mesh and evaluator (with symmetry support)
    model, symmetry_z0 = create_mesh_from_cfg(cfg, mesh_file)
    evaluator = build_evaluator_from_cfg(model, req, cfg; symmetry_z0=symmetry_z0)

    # Optimization config
    opt_config = build_opt_config(cfg)

    # Run
    result = run_metaheuristic(req, bounds, evaluator; config=opt_config)
    @info "Optimization completed. Best objective: $(result.best_value)"

    # Save best coils
    best_coils_json = joinpath(out_dir, "best_coils.json")
    coils_dict = [
        Dict(
            "r_inner" => c.r_inner,
            "r_outer" => c.r_outer,
            "height" => c.height,
            "z_center" => c.z_center,
            "current_density" => c.current_density
        ) for c in result.coils
    ]
    open(best_coils_json, "w") do f
        JSON.print(f, Dict("coils" => coils_dict, "objective" => result.best_value,
                          "base_mse" => result.base_mse,
                          "penalties" => Dict("field" => result.penalties.field,
                                             "geometry" => result.penalties.geometry,
                                             "safety" => result.penalties.safety,
                                             "conductor" => result.penalties.conductor),
                          "iterations" => result.iterations, "f_calls" => result.f_calls), 4)
    end

    # Field profile
    profile_json = joinpath(out_dir, "field_profile.json")
    z_samples = result.metrics.z_samples
    Bz_samples = result.metrics.Bz_samples
    profile_data = [
        Dict("z" => z, "Bz_computed" => Bz, "Bz_target" => req.target_profile.Bz[argmin(abs.(req.target_profile.z .- z))])
        for (z, Bz) in zip(z_samples, Bz_samples)
    ]
    open(profile_json, "w") do f
        JSON.print(f, Dict("profile" => profile_data), 4)
    end

    return (result=result, out_dir=out_dir)
end

end # module


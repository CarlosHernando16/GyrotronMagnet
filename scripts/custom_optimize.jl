"""
    custom_optimize.jl

Custom optimization script for GyrotronMagnet that allows defining arbitrary target field profiles
and constraints, then runs the coil geometry optimization.

Usage:
    julia --project=. scripts/custom_optimize.jl [--config CONFIG_FILE] [--output OUTPUT_DIR]

This script differs from optimize_geometry.jl by supporting custom target field profiles via JSON config.
Example configs are available in scripts/configs/.
It establishes requirements (target profile, geometric constraints, safety, penalties), sets up
FEM evaluation, and runs metaheuristic optimization.

Outputs: Same as optimize_geometry.jl (JSON files, VTK, console summary).
"""

using DrWatson
@quickactivate "GyrotronMagnet"

using ArgParse
using JSON
using Dates
using Gridap
using GyrotronMagnet

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config"
            help = "Path to JSON configuration file (optional, supports custom profiles)"
            arg_type = String
            default = nothing
        "--output"
            help = "Output directory (default: datadir(\"optimization\", \"custom_runs\", \"run_YYYY-MM-DD_HH-MM-SS\"))"
            arg_type = String
            default = nothing
        "--mesh-file"
            help = "Path to Gmsh mesh file (.msh)"
            arg_type = String
            default = nothing
        "--quick-test"
            help = "Run quick test with reduced iterations"
            action = :store_true
    end
    return parse_args(s)
end

function load_config_file(config_path::String)
    if !isfile(config_path)
        error("Configuration file not found: $config_path")
    end
    return JSON.parsefile(config_path)
end

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
            "r_inner" => [0.06, 0.15],     # Narrower radial bounds
            "r_outer" => [0.10, 0.25],     # Narrower radial bounds
            "height" => [0.02, 0.06],      # Narrower height bounds
            "z_center" => [-0.10, 0.40],   # Focus on z=0 to 0.5 region
            "current_density" => [2.0e7, 6.0e7]  # Realistic current density range
        ),
        "optimization" => Dict(
            "algorithm" => "de",
            "population" => 50,      # Larger population
            "iterations" => 50,      # Fewer iterations but better search
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

function config_to_requirements(cfg::Dict)
    # Custom target profile
    profile_dict = cfg["target_profile"]
    profile = TargetFieldProfile(Float64.(profile_dict["z"]), Float64.(profile_dict["Bz"]);
                                 symmetry=profile_dict["symmetry"])

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
        get(geom_dict, "max_current_density", 8.0e7)  # Default fallback
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

function create_mesh_from_config(cfg::Dict, mesh_file=nothing)
    if mesh_file !== nothing && isfile(mesh_file)
        @info "Loading mesh from file: $mesh_file"
        model, _ = load_mesh(mesh_file)
        return model
    end

    mesh_cfg = cfg["mesh"]
    domain = mesh_cfg["domain"]
    resolution = mesh_cfg["resolution"]

    @info "Creating structured mesh: domain=$domain, resolution=$resolution"
    r_range = Tuple(domain[1])
    z_range = Tuple(domain[2])
    n_r = max(10, Int(ceil((r_range[2] - r_range[1]) / resolution)))
    n_z = max(10, Int(ceil((z_range[2] - z_range[1]) / resolution)))

    # Convert nested tuple format to flat 4-tuple format expected by Gridap
    domain_flat = (r_range[1], r_range[2], z_range[1], z_range[2])
    model = CartesianDiscreteModel(domain_flat, (n_r, n_z))
    return model
end

function main()
    args = parse_commandline()

    # Load or create configuration
    if args["config"] !== nothing
        cfg = merge_with_defaults(load_config_file(args["config"]))
    else
        @info "No config file provided, using defaults"
        cfg = create_default_config()
    end

    # Adjust for quick test
    if args["quick-test"]
        cfg["optimization"]["iterations"] = 10
        cfg["optimization"]["population"] = 10
        @info "Quick test mode: reduced iterations and population"
    end

    mesh_file = args["mesh-file"]
    run_and_save(cfg; mesh_file=mesh_file, output_dir=args["output"], quick=args["quick-test"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

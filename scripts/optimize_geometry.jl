"""
    optimize_geometry.jl

Main workflow script for optimizing gyrotron magnet coil geometries.

Usage:
    julia optimize_geometry.jl [--config CONFIG_FILE] [--output OUTPUT_DIR]

The script performs:
1. Loads optimization requirements and parameter bounds
2. Sets up FEM evaluator with mesh and solver configuration
3. Runs metaheuristic optimization (DE by default)
4. Saves results to DrWatson-compliant output directory

Outputs:
- Best coil configuration (JSON)
- Optimization history (convergence data)
- Final field profile comparison
- VTK export of best solution
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
            help = "Path to JSON configuration file (optional)"
            arg_type = String
            default = nothing
        "--output"
            help = "Output directory (default: datadir(\"optimization\", \"runs\"))"
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
        "requirements" => Dict(
            "B_target_central" => 2.0,
            "homogeneity_region_z" => [0.0, 0.25],
            "homogeneity_tolerance_ppm" => 2500.0,
            "decay_target_z" => 0.195,
            "decay_target_percentage" => 0.165,
            "bore_diameter_warm" => 0.180,
            "magnet_height_max" => 0.200,
            "outer_diameter_max" => 0.500,
            "alignment_tolerance_z" => 0.0005
        ),
        "bounds" => Dict(
            "n_coils" => 4,
            "r_inner" => [0.05, 0.2],
            "r_outer" => [0.08, 0.35],
            "height" => [0.01, 0.08],
            "z_center" => [-0.25, 0.25],
            "current_density" => [2.0e7, 1.0e8]
        ),
        "optimization" => Dict(
            "algorithm" => "de",
            "population" => 30,
            "iterations" => 100,
            "seed" => 0x12345678,
            "verbose" => true
        ),
        "mesh" => Dict(
            "domain" => [[0.0, 0.5], [0.0, 1.0]],
            "resolution" => 0.02
        )
    )
end

function config_to_requirements(cfg::Dict)
    req_dict = cfg["requirements"]
    profile = default_target_profile()
    em = ElectromagneticRequirements(
        req_dict["B_target_central"],
        req_dict["homogeneity_tolerance_ppm"],
        Tuple(req_dict["homogeneity_region_z"]),
        (req_dict["decay_target_z"], req_dict["decay_target_percentage"])
    )
    geom = GeometryConstraints(
        req_dict["bore_diameter_warm"],
        (req_dict["magnet_height_max"] * 0.9, req_dict["magnet_height_max"]),
        req_dict["outer_diameter_max"],
        req_dict["alignment_tolerance_z"],
        (2, 12),
        get(req_dict, "max_current_density", 1.0e8)
    )
    safety = SafetyConstraints(8.0, 4.0, true)
    penalties = PenaltyWeights(1.0, 1.0, 1.0, 1.0)
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
    # Gridap expects: (r_min, r_max, z_min, z_max), not ((r_min, r_max), (z_min, z_max))
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
    
    run_and_save(cfg; mesh_file=args["mesh-file"], output_dir=args["output"], quick=args["quick-test"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


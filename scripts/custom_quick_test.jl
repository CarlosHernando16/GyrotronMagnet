using Revise
using DrWatson
@quickactivate "GyrotronMagnet"
using GyrotronMagnet

# Load all functions from the script
include("custom_optimize.jl")

# Create and run quick test manually using the shared pipeline
cfg = create_default_config()
cfg["optimization"]["iterations"] = 20  # Better quick test settings
cfg["optimization"]["population"] = 25

# Uncomment to load/merge a custom config file:
# cfg = merge_with_defaults(load_config_file("scripts/configs/config_simple.json"))

# --- customize bounds here ---
cfg["bounds"]["n_coils"] = [4, 6]                     # integer
cfg["bounds"]["r_inner"] = 0.12         # meters
cfg["bounds"]["r_outer"] = [0.15, 0.2]          # meters
cfg["bounds"]["height"]  = 0.01                 # meters
cfg["bounds"]["z_center"] = [0, 0.35]        # meters
cfg["bounds"]["current_density"] = [2.5e7, 5.5e7] # A/m²
cfg["mesh"]["symmetry_z0"] = true

# --- end customization ---

println("Starting quick optimization test...")
result = run_and_save(cfg; quick=true)
println("Quick test completed! Best objective: $(result.result.best_value)")
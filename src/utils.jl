using DrWatson
using CSV
using DataFrames
using Dates
using JSON
using Gridap: Point

"""
    extract_field_on_axis(B_z::CellField, z_points::Vector{Float64})

Extract B_z component along the z-axis (r=0) for axisymmetric 2D problems.

# Arguments
- `B_z::CellField`: Axial magnetic field component (scalar field)
- `z_points::Vector{Float64}`: Z coordinates where to extract field

# Returns
- `B_z_values::Vector{Float64}`: Axial field values at specified z points (tesla)

# Notes
- For axisymmetric problems, on-axis (r≈0) we extract B_z component
- Points are evaluated at r=1e-6 to avoid singularity at r=0
"""
function extract_field_on_axis(B_z::CellField, z_points::Vector{Float64})
    B_z_values = Float64[]
    
    for z in z_points
        # Point near axis: (r≈0, z); use Gridap Point for evaluation
        x = Point(1e-6, z)
        
        # Evaluate B_z field at point
        B_z_val = B_z(x)
        push!(B_z_values, B_z_val)
    end
    
    return B_z_values
end

"""
    compute_field_profile(B_r::CellField, B_z::CellField, B_mag::CellField, 
                         mesh::DiscreteModel, path::Vector{Tuple{Float64, Float64}})

Extract magnetic field components along an arbitrary path for axisymmetric 2D problems.

# Arguments
- `B_r::CellField`: Radial magnetic field component (scalar field)
- `B_z::CellField`: Axial magnetic field component (scalar field)
- `B_mag::CellField`: Field magnitude (scalar field)
- `mesh::DiscreteModel`: Gridap discrete model
- `path::Vector{Tuple{Float64, Float64}}`: Path as vector of (r, z) coordinates

# Returns
- `B_r_values::Vector{Float64}`: Radial field component along path (tesla)
- `B_z_values::Vector{Float64}`: Axial field component along path (tesla)
- `B_mag_values::Vector{Float64}`: Field magnitude along path (tesla)

# Example
```julia
path = [(0.0, 0.0), (0.0, 0.1), (0.0, 0.2)]  # Along z-axis
B_r, B_z, B_mag = compute_field_profile(B_r_field, B_z_field, B_mag_field, model, path)
```
"""
function compute_field_profile(B_r::CellField, B_z::CellField, B_mag::CellField,
                               mesh::DiscreteModel, path::Vector{Tuple{Float64, Float64}})
    B_r_values = Float64[]
    B_z_values = Float64[]
    B_mag_values = Float64[]
    
    for (r, z) in path
        # Use tuple for Gridap compatibility
        x = (r, z)
        push!(B_r_values, B_r(x))
        push!(B_z_values, B_z(x))
        push!(B_mag_values, B_mag(x))
    end
    
    return B_r_values, B_z_values, B_mag_values
end

"""
    save_results(A_φ::FEFunction, B_z::CellField, params::Dict, filename::String)

Save FEM solution results using DrWatson for axisymmetric 2D problems.

# Arguments
- `A_φ::FEFunction`: Scalar potential solution (azimuthal component)
- `B_z::CellField`: Axial magnetic field component
- `params::Dict`: Dictionary of parameters to save
- `filename::String`: Base filename (without extension)

# Notes
- Saves to `datadir()` using DrWatson conventions
- Saves as JSON (simplified - full version would use JLD2 for fields)
- Includes metadata from params dictionary
"""
function save_results(A_φ::FEFunction, B_z::CellField, params::Dict, filename::String)
    # Use DrWatson's datadir for data storage
    data_path = datadir(filename * ".json")
    mkpath(dirname(data_path))
    
    # Save parameters as JSON along with metadata
    # Full implementation could serialize FE fields via JLD2
    results = Dict(
        "parameters" => params,
        "timestamp" => string(now()),
        "solver" => "GyrotronMagnet Phase 1"
    )
    
    # Save to JSON (simplified - full version would use JLD2 for fields)
    open(data_path, "w") do f
        JSON.print(f, results, 4)
    end
    
    @info "Results saved to: $data_path"
    
    return data_path
end

"""
    load_target_profile(filename::String)

Load target magnetic field profile from CSV file.

# Arguments
- `filename::String`: Path to CSV file

# Returns
- `z::Vector{Float64}`: Z coordinates
- `B_target::Vector{Float64}`: Target B_z values

# CSV Format
Expected format: z, B_z
Header row optional.
"""
function load_target_profile(filename::String)
    if !isfile(filename)
        error("Target profile file not found: $filename")
    end
    
    # Load CSV
    df = CSV.read(filename, DataFrames.DataFrame)
    
    # Assume columns are z and B_z
    z = df[:, 1] |> Vector{Float64}
    B_target = df[:, 2] |> Vector{Float64}
    
    return z, B_target
end

"""
    compute_rms_error(B_computed::Vector{Float64}, B_target::Vector{Float64})

Compute RMS error between computed and target fields.

# Arguments
- `B_computed::Vector{Float64}`: Computed field values
- `B_target::Vector{Float64}`: Target field values

# Returns
- `rms_error::Float64`: Root mean square error

# Formula
RMS = sqrt(mean((B_computed - B_target)²))
"""
function compute_rms_error(B_computed::Vector{Float64}, B_target::Vector{Float64})
    if length(B_computed) != length(B_target)
        error("Computed and target vectors must have same length")
    end
    
    n = length(B_computed)
    squared_errors = [(B_computed[i] - B_target[i])^2 for i in 1:n]
    rms_error = sqrt(sum(squared_errors) / n)
    
    return rms_error
end

# Export functions
export extract_field_on_axis, compute_field_profile, save_results, 
       load_target_profile, compute_rms_error


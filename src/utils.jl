using DrWatson
using CSV
using DataFrames
using Dates
using JSON
using Gridap: Point, Triangulation, writevtk

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

"""
    FieldGrid

Cache of magnetic-field samples over a structured r-z grid.

# Fields
- `r_vals::Vector{Float64}`: Radial sampling coordinates (actual locations where fields were evaluated)
- `z_vals::Vector{Float64}`: Axial sampling coordinates
- `data::Dict{Symbol, Matrix{Float64}}`: Cached matrices keyed by `:magnitude`, `:B_r`, `:B_z`

# Notes
- Matrices are sized (`length(r_vals)`, `length(z_vals)`)
- Radial coordinates are clamped to ≥ 1e-6 to avoid on-axis singularities
- Reuse the same cache to avoid re-evaluating FEM fields when plotting multiple components
"""
struct FieldGrid
    r_vals::Vector{Float64}
    z_vals::Vector{Float64}
    data::Dict{Symbol, Matrix{Float64}}
end

# Additional plotting and export helpers
using Plots

function _extract_model_bounds(model::DiscreteModel)
    Ω = Triangulation(model)
    coords = Ω.grid.node_coords
    flat_coords = coords isa AbstractMatrix ? vec(coords) : collect(coords)
    
    if isempty(flat_coords)
        return (0.0, 1.0), (0.0, 1.0)
    end
    
    r_vals = Float64[p[1] for p in flat_coords]
    z_vals = Float64[p[2] for p in flat_coords]
    return (minimum(r_vals), maximum(r_vals)), (minimum(z_vals), maximum(z_vals))
end

"""
    sample_field_grid(solution::MagnetostaticModel; n_r=50, n_z=100, 
                      r_limits=nothing, z_limits=nothing)

Sample all magnetic-field components on a structured r-z grid.

# Keyword Arguments
- `n_r::Int`: Number of radial samples (default: 50)
- `n_z::Int`: Number of axial samples (default: 100)
- `r_limits::Union{Nothing,Tuple{Float64,Float64}}`: Optional (r_min, r_max)
- `z_limits::Union{Nothing,Tuple{Float64,Float64}}`: Optional (z_min, z_max)

# Returns
- `FieldGrid`: Cache with sampled data for magnitude, B_r, and B_z

# Notes
- The cache can be reused across multiple plotting calls to avoid recomputation
- Radial coordinates are clamped to ≥ 1e-6 to avoid on-axis singularities
- The returned FieldGrid.r_vals contains the actual coordinates where fields were sampled
"""
function sample_field_grid(solution::MagnetostaticModel; n_r::Int=50, n_z::Int=100,
                           r_limits=nothing, z_limits=nothing)
    default_r, default_z = _extract_model_bounds(solution.model)
    r_min, r_max = something(r_limits, default_r)
    z_min, z_max = something(z_limits, default_z)

    r_vals = collect(range(r_min, r_max; length=n_r))
    z_vals = collect(range(z_min, z_max; length=n_z))

    # Store the actual r coordinates where data was sampled (after clamping)
    r_sampled = [max(r, 1e-6) for r in r_vals]

    grids = Dict(
        :magnitude => zeros(n_r, n_z),
        :B_r => zeros(n_r, n_z),
        :B_z => zeros(n_r, n_z)
    )

    for (i, r_safe) in enumerate(r_sampled)
        for (j, z) in enumerate(z_vals)
            pt = Point(r_safe, z)
            grids[:B_r][i, j] = solution.B_r(pt)
            grids[:B_z][i, j] = solution.B_z(pt)
            grids[:magnitude][i, j] = solution.B(pt)
        end
    end

    return FieldGrid(r_sampled, z_vals, grids)
end

"""
    plot_mesh(model::DiscreteModel; title="Mesh Plot")

Generate a 2D plot of the mesh structure, drawing actual grid lines when available.

# Arguments
- `model::DiscreteModel`: Gridap discrete model
- `title::String`: Plot title (default: "Mesh Plot")
"""
function plot_mesh(model::DiscreteModel; title::String="Mesh Plot")
    Ω = Triangulation(model)
    coords = Ω.grid.node_coords
    p = plot(title=title, xlabel="r [m]", ylabel="z [m]", aspect_ratio=:equal, legend=false)
    
    if coords isa AbstractMatrix
        n_r, n_z = size(coords)
        for i in 1:n_r
            rs = [coords[i, j][1] for j in 1:n_z]
            zs = [coords[i, j][2] for j in 1:n_z]
            plot!(p, rs, zs, color=:black, linewidth=0.4, alpha=0.7)
        end
        for j in 1:n_z
            rs = [coords[i, j][1] for i in 1:n_r]
            zs = [coords[i, j][2] for i in 1:n_r]
            plot!(p, rs, zs, color=:black, linewidth=0.4, alpha=0.7)
        end
    else
        flat_coords = collect(coords)
        scatter!(p, Float64[pt[1] for pt in flat_coords], Float64[pt[2] for pt in flat_coords];
                 markersize=1.5, markercolor=:black, alpha=0.6)
    end
    return p
end

"""
    plot_b_field(solution::MagnetostaticModel; field_type=:magnitude, title="B Field Plot",
                 grid_data=nothing, n_r=50, n_z=100, r_limits=nothing, z_limits=nothing,
                 plot_kind=:contourf)

Generate a 2D plot of the requested magnetic-field component. Optionally reuse a
cached grid (`FieldGrid`) to speed up repeated plotting.
"""
function plot_b_field(solution::MagnetostaticModel; field_type::Symbol=:magnitude,
                      title::String="B Field Plot", grid_data::Union{Nothing,FieldGrid}=nothing,
                      n_r::Int=50, n_z::Int=100, r_limits=nothing, z_limits=nothing,
                      plot_kind::Symbol=:contourf)
    data_cache = isnothing(grid_data) ?
        sample_field_grid(solution; n_r=n_r, n_z=n_z, r_limits=r_limits, z_limits=z_limits) :
        grid_data
    
    values = get(data_cache.data, field_type) do
        available = join(string.(collect(keys(data_cache.data))), ", ")
        error("Field type $(field_type) not in cache. Available: $available")
    end
    
    label = field_type == :magnitude ? "B magnitude [T]" :
            field_type == :B_r ? "B_r [T]" :
            field_type == :B_z ? "B_z [T]" :
            error("Unsupported field_type $field_type")
    
    plotting_fn = plot_kind == :heatmap ? heatmap : contourf
    p = plotting_fn(data_cache.r_vals, data_cache.z_vals, values';
                    title=title,
                    xlabel="r [m]",
                    ylabel="z [m]",
                    colorbar_title=label,
                    aspect_ratio=:equal,
                    linewidth=plot_kind == :contourf ? 0 : nothing)
    return p
end

"""
    export_solution_vtk(solution::MagnetostaticModel, filename::String; output_dir=datadir())

Export FEM solution data to VTK (VTU) for external visualization tools (ParaView, VisIt, etc.).
"""
function export_solution_vtk(solution::MagnetostaticModel, filename::String; output_dir=datadir())
    mkpath(output_dir)
    basepath = joinpath(output_dir, filename)
    
    writevtk(solution.Ω, basepath;
             cellfields=[
                "A_phi" => solution.A_φ,
                "B_mag" => solution.B,
                "B_r" => solution.B_r,
                "B_z" => solution.B_z
             ])
    return basepath * ".vtu"
end

# Export functions
export extract_field_on_axis, compute_field_profile, save_results, 
       load_target_profile, compute_rms_error, FieldGrid, sample_field_grid,
       plot_mesh, plot_b_field, export_solution_vtk



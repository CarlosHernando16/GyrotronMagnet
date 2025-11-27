"""
    CoilParameters

Structure storing coil configuration parameters for axisymmetric 2D problems.

# Fields
- `r::Float64`: Radial position of coil center (meters)
- `z::Float64`: Axial position of coil center (meters)
- `I::Float64`: Current in coil (amps)
- `n_turns::Int`: Number of turns
- `r_inner::Float64`: Inner radius of coil in meters. These radii are absolute distances from the axis.
- `r_outer::Float64`: Outer radius of coil in meters (also absolute).
- `z_length::Float64`: Axial length of the coil (meters)

# Notes
- The coil occupies the region `r_inner ≤ r ≤ r_outer` and `z - z_length/2 ≤ z ≤ z + z_length/2`.
- Current density is computed as `J_φ = n_turns * I / (cross-section area)`.
"""
struct CoilParameters
    r::Float64      # Radial position of coil center
    z::Float64      # Axial position of coil center
    I::Float64      # Current (amps)
    n_turns::Int    # Number of turns
    r_inner::Float64
    r_outer::Float64
    z_length::Float64
end

"""
    create_coil_geometry(r, z, I, n_turns; r_inner=0.01, r_outer=0.02, z_length=0.01)

Create a coil geometry with specified parameters.

# Arguments
- `r::Float64`: Radial position of coil center in meters
- `z::Float64`: Axial position of coil center in meters
- `I::Float64`: Current in amps
- `n_turns::Int`: Number of turns
- `r_inner::Float64`: Inner radius in meters (default 0.01 m)
- `r_outer::Float64`: Outer radius in meters (default 0.02 m)
- `z_length::Float64`: Axial length in meters (default 0.01 m)

# Returns
- `CoilParameters`: Coil configuration structure
"""
function create_coil_geometry(r=0.1, z=0.5, I=1000.0, n_turns=100; r_inner=0.01, r_outer=0.02, z_length=0.01)
    return CoilParameters(r, z, I, n_turns, r_inner, r_outer, z_length)
end

"""
    generate_mesh(domain::Tuple, resolution::Float64; filename=nothing)

Generate a 2D axisymmetric mesh using GridapGmsh.

# Arguments
- `domain::Tuple`: Domain bounds as ((r_min, r_max), (z_min, z_max))
- `resolution::Float64`: Mesh resolution parameter
- `filename::Union{String, Nothing}`: Optional filename to save mesh (default: nothing)

# Returns
- `model::DiscreteModel`: Gridap discrete model
- `labels::FaceLabels`: Labels for boundary conditions

# Example
```julia
domain = ((0.0, 0.5), (0.0, 1.0))  # r: 0-0.5m, z: 0-1.0m
model, labels = generate_mesh(domain, 0.05)
```
"""
function generate_mesh(domain::Tuple, resolution::Float64; filename=nothing)
    (r_range, z_range) = domain
    r_min, r_max = r_range
    z_min, z_max = z_range
    
    # Create Gmsh model for 2D axisymmetric geometry
    # For now, create a simple rectangular domain
    # In practice, you would use Gmsh API to create more complex geometries
    
    # Generate mesh using GridapGmsh
    # This is a simplified version - full implementation would use Gmsh.jl
    # For now, we'll create a structured mesh approximation
    
    # Create a simple structured mesh
    # Note: GridapGmsh requires a .msh file, so we'll need to generate one
    # For Phase 1, we'll provide a basic implementation
    
    if filename === nothing
        filename = joinpath(projectdir(), "data", "mesh_temp.msh")
        mkpath(dirname(filename))
    end
    
    # Create Gmsh geometry script
    gmsh_script = """
    // 2D axisymmetric mesh for gyrotron magnet
    SetFactory("OpenCASCADE");
    
    // Domain rectangle
    Rectangle(1) = {$r_min, $z_min, 0, $(r_max-r_min), $(z_max-z_min), 0};
    
    // Mesh resolution
    Mesh.CharacteristicLengthMin = $resolution;
    Mesh.CharacteristicLengthMax = $resolution;
    
    // Generate mesh
    Mesh 2;
    
    // Save mesh
    Save "$filename";
    """
    
    # Write script to temporary file
    script_file = filename * ".geo"
    open(script_file, "w") do f
        write(f, gmsh_script)
    end
    
    # Note: This requires Gmsh to be installed and available
    # For Phase 1, we'll use GridapGmsh to load the mesh
    # In practice, you would run: gmsh script_file -2 -o filename
    
    # Load mesh using GridapGmsh
    # For now, return a placeholder that indicates mesh needs to be generated
    # Full implementation requires Gmsh.jl or manual mesh file creation
    
    error("Mesh generation requires Gmsh installation. Please use load_mesh() with an existing .msh file, or install Gmsh and implement full mesh generation.")
end

"""
    load_mesh(msh_file::String)

Load an existing Gmsh mesh file.

# Arguments
- `msh_file::String`: Path to .msh file

# Returns
- `model::DiscreteModel`: Gridap discrete model
- `labels::FaceLabels`: Labels for boundary conditions

# Example
```julia
model, labels = load_mesh("data/mesh.msh")
```
"""
function load_mesh(msh_file::String)
    if !isfile(msh_file)
        error("Mesh file not found: $msh_file")
    end
    
    # Load mesh using GridapGmsh
    model = GmshDiscreteModel(msh_file)
    
    # Get labels for boundary conditions
    labels = get_face_labeling(model)
    
    return model, labels
end

"""
    create_current_density(coil_params::Vector{CoilParameters}, model::DiscreteModel)

Create azimuthal current density field J_φ from coil parameters for 2D axisymmetric problems.

# Arguments
- `coil_params::Vector{CoilParameters}`: Vector of coil configurations
- `model::DiscreteModel`: Gridap discrete model (2D r-z plane)

# Returns
- `J_φ::CellField`: Azimuthal current density (scalar field) in A/m²

# Notes
- For axisymmetric 2D, only J_φ (azimuthal) component exists
- Current density is computed as: J_φ = n_turns * I / (coil cross-section area)
- The area is: π * (r_outer² - r_inner²) * z_length
- Returns a scalar field (not vector) for axisymmetric formulation
"""
function create_current_density(coil_params::Vector{CoilParameters}, model::DiscreteModel)
    # Get domain
    Ω = Triangulation(model)
    
    # Create current density function (scalar J_φ)
    function J_φ_field(x)
        r = x[1]  # Radial coordinate
        z = x[2]  # Axial coordinate
        
        J_total = 0.0
        
        # Sum contributions from all coils
        for coil in coil_params
            # Check if point is inside coil
            # For axisymmetric, coil is an annulus: r_inner <= r <= r_outer
            # and z_cent - z_length/2 <= z <= z_cent + z_length/2
            z_low = coil.z - coil.z_length/2
            z_high = coil.z + coil.z_length/2
            
            if (coil.r_inner <= r <= coil.r_outer) && (z_low <= z <= z_high)
                # Compute current density: J_φ = n_turns * I / area
                # Area of annulus: π * (r_outer² - r_inner²) * z_length
                coil_area = π * (coil.r_outer^2 - coil.r_inner^2) * coil.z_length
                J_total += coil.n_turns * coil.I / coil_area
            end
        end
        
        return J_total
    end
    
    # Create CellField from function (scalar field)
    J_φ = CellField(J_φ_field, Ω)
    
    return J_φ
end

# Export functions and types
export CoilParameters, create_coil_geometry, generate_mesh, load_mesh, create_current_density

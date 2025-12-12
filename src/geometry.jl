"""
    CoilGeometry

Geometry and source term parameters for a single axisymmetric coil.

# Fields
- `r::Float64`: Radial position of coil center (meters)
- `z::Float64`: Axial position of coil center (meters)
- `r_inner::Float64`: Inner radius of coil in meters (absolute distance from axis)
- `r_outer::Float64`: Outer radius of coil in meters (absolute distance from axis)
- `z_length::Float64`: Axial length of the coil (meters)
- `current_density::Float64`: Azimuthal current density J_φ in A/m² (FEM source term)
"""
struct CoilGeometry
    r::Float64
    z::Float64
    r_inner::Float64
    r_outer::Float64
    z_length::Float64
    current_density::Float64
end

"""
    ConductorParameters

Engineering properties of the conductor used to wind a coil.

# Fields
- `wire_cross_section::Float64`: Wire cross-sectional area (m²)
- `material::String`: Conductor material name
- `winding_factor::Float64`: Packing efficiency (0-1), accounts for voids/insulation
- `insulation_thickness::Float64`: Electrical insulation thickness (m)
"""
struct ConductorParameters
    wire_cross_section::Float64
    material::String
    winding_factor::Float64
    insulation_thickness::Float64
    function ConductorParameters(wire_cross_section::Real,
                                 material::AbstractString;
                                 winding_factor::Real = 0.7,
                                 insulation_thickness::Real = 0.0)
        wc = float(wire_cross_section)
        wf = float(winding_factor)
        it = float(insulation_thickness)
        wc > 0 || error("wire_cross_section must be positive.")
        (0 < wf <= 1) || error("winding_factor must be in (0, 1].")
        it >= 0 || error("insulation_thickness must be non-negative.")
        return new(wc, String(material), wf, it)
    end
end

"""
    create_coil_geometry(; r, z, current_density, r_inner=0.01, r_outer=0.02, z_length=0.01)

Create a coil geometry with specified parameters using current density as the source term.
All arguments are keywords for consistency.

# Arguments
- `r::Real`: Radial position of coil center in meters
- `z::Real`: Axial position of coil center in meters
- `current_density::Real`: Azimuthal current density J_φ in A/m²
- `r_inner::Real`: Inner radius in meters (default 0.01 m)
- `r_outer::Real`: Outer radius in meters (default 0.02 m)
- `z_length::Real`: Axial length in meters (default 0.01 m)

# Returns
- `CoilGeometry`: Coil configuration structure
"""
function create_coil_geometry(; r::Real, z::Real, current_density::Real,
                               r_inner::Real=0.01, r_outer::Real=0.02, z_length::Real=0.01)
    return CoilGeometry(float(r), float(z), float(r_inner), float(r_outer),
                        float(z_length), float(current_density))
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
    create_current_density(coil_geometries::Vector{CoilGeometry}, model::DiscreteModel)

Create azimuthal current density field J_φ from coil parameters for 2D axisymmetric problems.

# Arguments
- `coil_geometries::Vector{CoilGeometry}`: Vector of coil configurations
- `model::DiscreteModel`: Gridap discrete model (2D r-z plane)

# Returns
- `J_φ::CellField`: Azimuthal current density (scalar field) in A/m²

# Notes
- For axisymmetric 2D, only J_φ (azimuthal) component exists.
- Current density is provided directly by `CoilGeometry.current_density` (already in A/m²);
  no turns or transport current are inferred here.
- Returns a scalar field (not vector) for axisymmetric formulation.
"""
function create_current_density(coil_geometries::Vector{CoilGeometry}, model::DiscreteModel)
    # Get domain
    Ω = Triangulation(model)
    
    # Create current density function (scalar J_φ)
    function J_φ_field(x)
        r = x[1]  # Radial coordinate
        z = x[2]  # Axial coordinate
        
        J_total = 0.0
        
        # Sum contributions from all coils
        for coil in coil_geometries
            # Check if point is inside coil
            # For axisymmetric, coil is an annulus: r_inner <= r <= r_outer
            # and z_cent - z_length/2 <= z <= z_cent + z_length/2
            z_low = coil.z - coil.z_length/2
            z_high = coil.z + coil.z_length/2
            
            if (coil.r_inner <= r <= coil.r_outer) && (z_low <= z <= z_high)
                # Current density is provided directly as design variable
                J_total += coil.current_density
            end
        end
        
        return J_total
    end
    
    # Create CellField from function (scalar field)
    J_φ = CellField(J_φ_field, Ω)
    
    return J_φ
end

"""
    compute_coil_area(geom::CoilGeometry)

Cross-sectional area of the coil pack (annulus) times axial length.
"""
compute_coil_area(geom::CoilGeometry) = π * (geom.r_outer^2 - geom.r_inner^2) * geom.z_length

"""
    compute_mean_circumference(geom::CoilGeometry)

Mean circumference of the annulus where turns are wound.
"""
compute_mean_circumference(geom::CoilGeometry) = π * (geom.r_inner + geom.r_outer)

"""
    compute_ampere_turns(geom::CoilGeometry)

Total ampere-turns implied by current density over the coil cross-section.
"""
compute_ampere_turns(geom::CoilGeometry) = geom.current_density * compute_coil_area(geom)

"""
    compute_total_current(geom::CoilGeometry)

Total current if interpreted as single-turn equivalent (A) from J_φ * area.
"""
compute_total_current(geom::CoilGeometry) = compute_ampere_turns(geom)

"""
    compute_n_turns(geom::CoilGeometry; operating_current)

Compute number of turns assuming a chosen operating current per turn.
"""
function compute_n_turns(geom::CoilGeometry; operating_current::Real)
    I = float(operating_current)
    I > 0 || error("operating_current must be positive.")
    return compute_ampere_turns(geom) / I
end

"""
    compute_conductor_length(geom::CoilGeometry, conductor::ConductorParameters; operating_current)

Estimate total conductor length using mean circumference, winding factor, and derived turns.
"""
function compute_conductor_length(geom::CoilGeometry, conductor::ConductorParameters; operating_current::Real)
    n_turns = compute_n_turns(geom; operating_current=operating_current)
    circ = compute_mean_circumference(geom)
    return (n_turns * circ) / conductor.winding_factor
end

# Export functions and types
export CoilGeometry, ConductorParameters, create_coil_geometry, generate_mesh, load_mesh,
       create_current_density, compute_coil_area, compute_mean_circumference,
       compute_ampere_turns, compute_total_current, compute_n_turns, compute_conductor_length

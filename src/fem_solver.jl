"""
    MagnetostaticModel

Structure storing magnetostatic FEM model and solution for 2D axisymmetric problems.

# Fields
- `model::DiscreteModel`: Gridap discrete model
- `Ω::Triangulation`: Domain triangulation
- `V::FESpace`: Finite element space for scalar potential A_φ
- `A_φ::FEFunction`: Solution scalar potential (azimuthal component)
- `B::CellField`: Computed magnetic field magnitude
- `B_r::CellField`: Radial component of magnetic field
- `B_z::CellField`: Axial component of magnetic field
"""
struct MagnetostaticModel
    model::DiscreteModel
    Ω::Triangulation
    V::FESpace
    A_φ::FEFunction
    B::CellField
    B_r::CellField
    B_z::CellField
end

"""
    solve_magnetostatic(model::DiscreteModel, J_φ::CellField;
                       μ₀::Float64=4π*1e-7, order::Int=1, bc_type="dirichlet",
                       boundary_tag="outer_boundary", symmetry_z0::Bool=false,
                       coil_geometries::Vector=[], mu_r_air::Float64=1.0, mu_r_coil::Float64=1.0)

Solve the 2D axisymmetric magnetostatic problem using scalar potential A_φ.

The governing equation in axisymmetric coordinates is:
∇×(ν∇×A_φ) = J_φ

which reduces to the scalar equation:
-(1/r)∇·(r/μ ∇A_φ) = J_φ

# Arguments
- `model::DiscreteModel`: Gridap discrete model (2D r-z plane)
- `J_φ::CellField`: Azimuthal current density in A/m²
- `μ₀::Float64`: Permeability of free space (default 4π×10⁻⁷ H/m)
- `order::Int`: FE polynomial order (default 1)
- `bc_type::String`: Boundary condition type (default `"dirichlet"`)
- `boundary_tag::String`: Tag name for boundary conditions (default `"outer_boundary"`)
- `symmetry_z0::Bool`: If true, apply Neumann BC at z=0 (symmetry plane) and Dirichlet at outer boundaries
- `coil_geometries::Vector`: Vector of `CoilGeometry` for region-aware material properties
- `mu_r_air::Float64`: Relative permeability of air region (default 1.0)
- `mu_r_coil::Float64`: Relative permeability of coil region (default 1.0)

# Returns
- `MagnetostaticModel`: Model containing solution A_φ and computed B fields

# Mathematical Formulation
The weak form for axisymmetric 2D is:
∫_Ω (ν(x)) * (r * ∇u · ∇v + (u*v)/r) dΩ = ∫_Ω J_φ * v * r dΩ, ∀v ∈ V₀

where:
- u = A_φ (azimuthal component of vector potential)
- r is the radial coordinate
- ν(x) = 1/(μ₀ * μᵣ(x)) is the reluctivity (spatially varying if materials differ)
- B_r = -∂A_φ/∂z
- B_z = (1/r)∂(r*A_φ)/∂r

# Symmetry Mode
- r = 0 (axis): Always Dirichlet (A_φ = 0) for axisymmetric regularity.
- Outer boundaries (r_max, z_min, z_max): Dirichlet (A_φ = 0).
- When `symmetry_z0=true`, we solve only z ≥ 0; z=0 uses Neumann (∂A/∂n = 0) to enforce
  midplane symmetry, while keeping r=0 Dirichlet.

# Notes
- Domain is in r-z plane (axisymmetric)
- Only A_φ component is non-zero
- Magnetic field components: B_r (radial) and B_z (axial)
"""
function solve_magnetostatic(model::DiscreteModel, J_φ::CellField;
                             μ₀::Float64=4π*1e-7, order::Int=1, bc_type="dirichlet",
                             boundary_tag="outer_boundary", symmetry_z0::Bool=false,
                             coil_geometries::Vector=[], mu_r_air::Float64=1.0, mu_r_coil::Float64=1.0)
    
    # Get domain triangulation
    Ω = Triangulation(model)
    
    # Define FE space for scalar potential A_φ
    reffe = ReferenceFE(lagrangian, Float64, order)
    
    # Get face labeling for boundary conditions
    labels = get_face_labeling(model)
    tag_names = collect(keys(labels.tag_to_entities))
    
    # Apply boundary conditions
    if bc_type == "dirichlet"
        if symmetry_z0
            # Half-domain symmetry: Dirichlet on r=0, r_max, z_max; Neumann on z=0
            if boundary_tag in tag_names
                V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags=boundary_tag)
                U = TrialFESpace(V, 0.0)
            else
                V = TestFESpace(model, reffe, conformity=:H1,
                    dirichlet_masks=[(true, true, false, true)])  # (r_min, r_max, z_min, z_max)
                U = TrialFESpace(V, 0.0)
            end
        else
            # Standard full-domain: Dirichlet on all boundaries (including axis)
            if boundary_tag in tag_names
                V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags=boundary_tag)
                U = TrialFESpace(V, 0.0)
            else
                V = TestFESpace(model, reffe, conformity=:H1,
                    dirichlet_masks=[(true, true, true, true)])  # (r_min, r_max, z_min, z_max)
                U = TrialFESpace(V, 0.0)
            end
        end
    else
        # Pure Neumann (no Dirichlet constraints)
        V = TestFESpace(model, reffe, conformity=:H1)
        U = TrialFESpace(V)
    end
    
    # Define measure with appropriate quadrature order
    dΩ = Measure(Ω, 2*order)
    
    # Build spatially-varying reluctivity ν(x) = 1/(μ₀ * μᵣ(x))
    # If mu_r_air == mu_r_coil, use uniform reluctivity for efficiency
    if isempty(coil_geometries) || mu_r_air == mu_r_coil
        # Uniform material properties
        ν = 1.0 / (μ₀ * mu_r_air)
        ν_field = CellField(x -> ν, Ω)
    else
        # Spatially varying: check if point is inside any coil
        function compute_nu(x)
            r, z = x[1], x[2]
            for coil in coil_geometries
                z_low = coil.z - coil.z_length / 2
                z_high = coil.z + coil.z_length / 2
                if (coil.r_inner <= r <= coil.r_outer) && (z_low <= z <= z_high)
                    return 1.0 / (μ₀ * mu_r_coil)
                end
            end
            return 1.0 / (μ₀ * mu_r_air)
        end
        ν_field = CellField(compute_nu, Ω)
    end
    
    # Weak form for axisymmetric 2D magnetostatics with variable reluctivity
    # ∫ ν(x) * (r * ∇u · ∇v + (u*v)/r) dΩ = ∫ J_φ * v * r dΩ
    function a(u, v)
        r_coord = CellField(x -> x[1], Ω)
        return ∫(ν_field * (r_coord * (∇(u) ⋅ ∇(v)) + (u*v)/r_coord))dΩ
    end
    
    function l(v)
        r_coord = CellField(x -> x[1], Ω)
        return ∫(J_φ * v * r_coord)dΩ
    end
    
    # Assemble and solve
    op = AffineFEOperator(a, l, U, V)
    A_φ_sol = solve(op)
    
    # Compute magnetic field components
    # B_r = -∂A_φ/∂z
    # B_z = (1/r)∂(r*A_φ)/∂r = ∂A_φ/∂r + A_φ/r
    dA_φ = ∇(A_φ_sol)
    
    # B_r = -∂A_φ/∂z (negative of z-component of gradient)
    B_r_field = CellField(x -> -dA_φ(x)[2], Ω)
    
    # B_z = ∂A_φ/∂r + A_φ/r = dA_φ[1] + A_φ/r
    function B_z_func(x)
        r = x[1]
        if r < 1e-10  # On axis, use limit: B_z = 2*∂A_φ/∂r
            return 2 * dA_φ(x)[1]
        else
            return dA_φ(x)[1] + A_φ_sol(x) / r
        end
    end
    B_z_field = CellField(B_z_func, Ω)
    
    # Total field magnitude
    B_mag_field = CellField(x -> sqrt(B_r_field(x)^2 + B_z_field(x)^2), Ω)
    
    # Create and return model structure
    return MagnetostaticModel(model, Ω, V, A_φ_sol, B_mag_field, B_r_field, B_z_field)
end

# Backward-compatible method without keyword arguments for μ₀
function solve_magnetostatic(model::DiscreteModel, J_φ::CellField, μ₀::Float64; kwargs...)
    return solve_magnetostatic(model, J_φ; μ₀=μ₀, kwargs...)
end

"""
    compute_magnetic_field(A_φ::FEFunction, model::DiscreteModel)

Compute magnetic field components from scalar potential A_φ.

# Arguments
- `A_φ::FEFunction`: Scalar potential solution (azimuthal component)
- `model::DiscreteModel`: Gridap discrete model

# Returns
- `B_r::CellField`: Radial component of magnetic field
- `B_z::CellField`: Axial component of magnetic field
- `B_mag::CellField`: Field magnitude

# Notes
- B_r = -∂A_φ/∂z
- B_z = ∂A_φ/∂r + A_φ/r
- On axis (r=0): B_z = 2*∂A_φ/∂r
"""
function compute_magnetic_field(A_φ::FEFunction, model::DiscreteModel)
    Ω = Triangulation(model)
    dA_φ = ∇(A_φ)
    
    # B_r = -∂A_φ/∂z
    B_r = CellField(x -> -dA_φ(x)[2], Ω)
    
    # B_z = ∂A_φ/∂r + A_φ/r
    function B_z_func(x)
        r = x[1]
        if r < 1e-10
            return 2 * dA_φ(x)[1]
        else
            return dA_φ(x)[1] + A_φ(x) / r
        end
    end
    B_z = CellField(B_z_func, Ω)
    
    # Magnitude
    B_mag = CellField(x -> sqrt(B_r(x)^2 + B_z(x)^2), Ω)
    
    return B_r, B_z, B_mag
end

"""
    apply_boundary_conditions(V::FESpace, model::DiscreteModel, 
                             bc_type::String="dirichlet", boundary_tag="outer_boundary")

Apply boundary conditions to the FE space (helper function, usually called internally).

# Arguments
- `V::FESpace`: Finite element space
- `model::DiscreteModel`: Gridap discrete model
- `bc_type::String`: Type of boundary condition ("dirichlet" or "neumann")
- `boundary_tag::String`: Tag name for boundaries

# Returns
- `V0::FESpace`: Constrained FE space with boundary conditions applied

# Notes
- Dirichlet BC: A_φ = 0 on boundaries
- For axisymmetric problems, this is the standard boundary condition
"""
function apply_boundary_conditions(V::FESpace, model::DiscreteModel, 
                                   bc_type::String="dirichlet", boundary_tag="outer_boundary")
    if bc_type == "dirichlet"
        # Apply homogeneous Dirichlet BC: A_φ = 0 on boundaries
        V0 = TestFESpace(model, V.reffe, conformity=:H1, dirichlet_tags=boundary_tag)
        return V0
    else
        # Neumann BC: no constraints needed
        return V
    end
end

# Export functions and types
export MagnetostaticModel, solve_magnetostatic, compute_magnetic_field, apply_boundary_conditions

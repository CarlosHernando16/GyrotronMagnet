module OptimizationSpecs

"""
    TargetFieldProfile(z, Bz; symmetry=true)

Container for the desired on-axis magnetic field profile.

# Arguments
- `z::AbstractVector{<:Real}`: Axial sampling positions (meters)
- `Bz::AbstractVector{<:Real}`: Target axial field values (tesla)
- `symmetry::Bool`: Whether the profile must be symmetric around the equatorial plane (default `true`)

# Notes
- Vectors must be the same length and strictly sorted in `z`
- Values are promoted to `Float64` for stable downstream computations
"""
struct TargetFieldProfile
    z::Vector{Float64}
    Bz::Vector{Float64}
    symmetry::Bool
    function TargetFieldProfile(z::AbstractVector{<:Real}, Bz::AbstractVector{<:Real}; symmetry::Bool=true)
        length(z) == length(Bz) || error("Target profile z and Bz must have the same length.")
        length(z) >= 2 || error("Target profile must contain at least two samples.")
        issorted(z) || error("Target profile z samples must be sorted.")
        return new(collect(float.(z)), collect(float.(Bz)), symmetry)
    end
end

"""
    ElectromagneticRequirements(; center_field=2.0, homogeneity_ppm=2500,
                               homogeneity_window=(0.15, 0.25),
                               decay_point=(0.195, 0.16))

Design-time electromagnetic constraints.

# Keywords
- `center_field`: Required `B_z` at the cavity center (tesla)
- `homogeneity_ppm`: Max deviation within the homogeneous region (ppm)
- `homogeneity_window`: `(z_start, z_stop)` window (meters)
- `decay_point`: `(z_offset, fraction)` specifying desired decay at a given offset
"""
struct ElectromagneticRequirements
    center_field::Float64
    homogeneity_ppm::Float64
    homogeneity_window::Tuple{Float64,Float64}
    decay_point::Tuple{Float64,Float64}
end

"""
    GeometryConstraints(; bore_diameter=0.18, magnet_height=(0.18, 0.2),
                        outer_diameter=0.6, alignment_tol=5e-4,
                        coil_count=(2, 12), max_current_density=8.0e7)

Geometric and mechanical envelope limits.

# Keywords
- `bore_diameter`: Minimum clear bore (meters)
- `magnet_height`: `(min, max)` total height bounds (meters)
- `outer_diameter`: Maximum allowed outer diameter (meters)
- `alignment_tol`: Allowed axial misalignment from cavity (meters)
- `coil_count`: `(min, max)` number of coil groups
- `max_current_density`: Maximum conductor current density (A/m²)
"""
struct GeometryConstraints
    bore_diameter::Float64
    magnet_height::Tuple{Float64,Float64}
    outer_diameter::Float64
    alignment_tol::Float64
    coil_count::Tuple{Int,Int}
    max_current_density::Float64
end

"""
    SafetyConstraints(; max_conductor_field=8.0,
                       max_radial_field=4.0,
                       forbid_ferromagnetics=true)

Field and materials limits to preserve conductor safety.
"""
struct SafetyConstraints
    max_conductor_field::Float64
    max_radial_field::Float64
    forbid_ferromagnetics::Bool
end

"""
    PenaltyWeights(; field=1.0, geometry=1.0, safety=1.0, conductor=1.0)

Multipliers used to weight penalty terms in the optimization objective.
"""
struct PenaltyWeights
    field::Float64
    geometry::Float64
    safety::Float64
    conductor::Float64
end

"""
    OptimizationRequirements(profile, em_req, geom_req, safety, penalties)

Aggregate structure with all optimization inputs.
"""
struct OptimizationRequirements
    target_profile::TargetFieldProfile
    em::ElectromagneticRequirements
    geometry::GeometryConstraints
    safety::SafetyConstraints
    penalties::PenaltyWeights
end

"""
    default_target_profile(; z_max=0.5, dz=0.01)

Generate a symmetric default profile using the Phase 2 specs.
"""
function default_target_profile(; z_max=0.5, dz=0.01)
    z = collect(0.0:dz:z_max)
    Bz = similar(z)
    for (i, zz) in pairs(z)
        if zz ≤ 0.15
            Bz[i] = 2.0
        elseif zz ≤ 0.25
            frac = (zz - 0.15) / 0.10
            Bz[i] = 2.0 * (1 - 0.0025 * frac)
        else
            Bz[i] = 2.0 * 0.17^(zz / 0.195)
        end
    end
    return TargetFieldProfile(vcat(reverse(-z[2:end]), z), vcat(reverse(Bz[2:end]), Bz))
end

"""
    default_requirements()

Return an `OptimizationRequirements` instance populated with Phase 2 defaults.
"""
function default_requirements()
    profile = default_target_profile()
    em = ElectromagneticRequirements(2.0, 2500, (0.15, 0.25), (0.195, 0.17))
    geometry = GeometryConstraints(0.18, (0.18, 0.20), 0.60, 5e-4, (2, 12), 8.0e7)
    safety = SafetyConstraints(8.0, 4.0, true)
    penalties = PenaltyWeights(1.0, 1.0, 1.0, 1.0)
    return OptimizationRequirements(profile, em, geometry, safety, penalties)
end

"""
    validate_requirements(req::OptimizationRequirements)

Validate consistency of all requirements. Returns `true` on success, otherwise throws an `ArgumentError`.
"""
function validate_requirements(req::OptimizationRequirements)
    req.geometry.bore_diameter > 0 || error("Bore diameter must be positive.")
    req.geometry.magnet_height[1] < req.geometry.magnet_height[2] ||
        error("Magnet height bounds are invalid.")
    req.em.homogeneity_window[1] < req.em.homogeneity_window[2] ||
        error("Homogeneity window must be ordered.")
    req.geometry.coil_count[1] ≥ 1 || error("Minimum coil count must be ≥ 1.")
    req.geometry.max_current_density > 0 || error("Max current density must be positive.")
    req.safety.max_conductor_field > 0 || error("Safety field limits must be positive.")
    return true
end

# Exports
export TargetFieldProfile, ElectromagneticRequirements, GeometryConstraints,
       SafetyConstraints, PenaltyWeights, OptimizationRequirements,
       default_target_profile, default_requirements, validate_requirements

end # module


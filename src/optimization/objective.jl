module OptimizationObjective

using ..OptimizationSpecs: TargetFieldProfile, ElectromagneticRequirements,
    GeometryConstraints, SafetyConstraints, PenaltyWeights, OptimizationRequirements

"""
    CoilDescriptor(r_inner, r_outer, height, z_center, current_density)

Describes a single axisymmetric coil pack used during optimization.
"""
struct CoilDescriptor
    r_inner::Float64
    r_outer::Float64
    height::Float64
    z_center::Float64
    current_density::Float64  # Current density in A/m² (design variable)
end

"""
    SimulationMetrics(z_samples, Bz_samples; max_conductor_field, max_radial_field)

Stores the relevant simulation outputs required by the objective function.
"""
struct SimulationMetrics
    z_samples::Vector{Float64}
    Bz_samples::Vector{Float64}
    max_conductor_field::Float64
    max_radial_field::Float64
    function SimulationMetrics(z_samples::AbstractVector{<:Real},
                               Bz_samples::AbstractVector{<:Real};
                               max_conductor_field::Real,
                               max_radial_field::Real)
        length(z_samples) == length(Bz_samples) ||
            error("z_samples and Bz_samples must have equal length.")
        length(z_samples) >= 2 || error("Simulation metrics require at least two samples.")
        issorted(z_samples) || error("Simulation z_samples must be sorted.")
        return new(collect(float.(z_samples)),
                   collect(float.(Bz_samples)),
                   float(max_conductor_field),
                   float(max_radial_field))
    end
end

"""
    PenaltyBreakdown(field, geometry, safety, conductor)

Convenience type returning each penalty component.
"""
struct PenaltyBreakdown
    field::Float64
    geometry::Float64
    safety::Float64
    conductor::Float64
end

# --- helper utilities -------------------------------------------------------

function _resample(z_src::Vector{Float64}, v_src::Vector{Float64}, z_tgt::Vector{Float64})
    res = similar(z_tgt)
    for (i, z) in pairs(z_tgt)
        if z <= z_src[1]
            res[i] = v_src[1]
        elseif z >= z_src[end]
            res[i] = v_src[end]
        else
            j = searchsortedlast(z_src, z)
            t = (z - z_src[j]) / (z_src[j+1] - z_src[j])
            res[i] = (1 - t) * v_src[j] + t * v_src[j+1]
        end
    end
    return res
end

function _field_error(profile::TargetFieldProfile, metrics::SimulationMetrics)
    sim_on_target = _resample(metrics.z_samples, metrics.Bz_samples, profile.z)
    diff = sim_on_target .- profile.Bz
    return sum(abs2, diff), sim_on_target
end

function _homogeneity_penalty(profile_vals::Vector{Float64},
                              profile::TargetFieldProfile,
                              em::ElectromagneticRequirements)
    window = em.homogeneity_window
    idx = findall(z -> window[1] ≤ z ≤ window[2], profile.z)
    isempty(idx) && return 0.0
    # Guard against division by zero
    abs(em.center_field) < 1e-15 && return 0.0
    slice = profile_vals[idx]
    deviations = abs.((slice .- em.center_field) / em.center_field * 1e6)
    exceed = maximum(deviations) - em.homogeneity_ppm
    return exceed > 0 ? exceed^2 : 0.0
end

function _decay_penalty(metrics::SimulationMetrics, em::ElectromagneticRequirements)
    pt = _resample(metrics.z_samples, metrics.Bz_samples, [em.decay_point[1]])[1]
    target = em.center_field * em.decay_point[2]
    return abs(pt - target)^2
end

function _geometry_violations(coils::Vector{CoilDescriptor}, geom::GeometryConstraints)
    isempty(coils) && return (0.0, 0.0, 0.0, 0.0)
    min_inner = minimum(c.r_inner for c in coils)
    max_outer = maximum(c.r_outer for c in coils)
    z_low = minimum(c.z_center - c.height/2 for c in coils)
    z_high = maximum(c.z_center + c.height/2 for c in coils)
    total_height = z_high - z_low
    max_j = maximum(c.current_density for c in coils)
    v_bore = max(geom.bore_diameter/2 - min_inner, 0)
    v_outer = max(max_outer - geom.outer_diameter/2, 0)
    low, high = geom.magnet_height
    v_height = total_height < low ? low - total_height :
               total_height > high ? total_height - high : 0.0
    v_j = max(max_j - geom.max_current_density, 0)
    return (v_bore, v_outer, v_height, v_j)
end

function _safety_penalty(metrics::SimulationMetrics, safety::SafetyConstraints)
    p1 = max(metrics.max_conductor_field - safety.max_conductor_field, 0)
    p2 = max(metrics.max_radial_field - safety.max_radial_field, 0)
    return p1^2 + p2^2
end

function compute_penalties(profile_vals::Vector{Float64},
                           profile::TargetFieldProfile,
                           metrics::SimulationMetrics,
                           req::OptimizationRequirements,
                           coils::Vector{CoilDescriptor})
    field_pen = _homogeneity_penalty(profile_vals, profile, req.em) +
                _decay_penalty(metrics, req.em)
    
    bore_v, outer_v, height_v, j_v = _geometry_violations(coils, req.geometry)
    geom_pen = bore_v^2 + outer_v^2 + height_v^2 + j_v^2
    
    safety_pen = _safety_penalty(metrics, req.safety)

    # Note: Conductor current density penalty removed because:
    # - c.current_density is FEM azimuthal current density (A/m²), not conductor wire current density
    # - Conductor safety is already constrained by magnetic field limits in safety_pen
    # - Wire cross-sectional information is not available in CoilDescriptor
    conductor_pen = 0.0

    return PenaltyBreakdown(field_pen, geom_pen, safety_pen, conductor_pen)
end

"""
    objective_function(metrics, coils, req)

Compute the total cost (field mismatch + penalties) given simulation metrics and coil descriptors.

Returns `(total, base_mse, penalties::PenaltyBreakdown)`.
"""
function objective_function(metrics::SimulationMetrics,
                            coils::Vector{CoilDescriptor},
                            req::OptimizationRequirements)
    base_mse, sim_on_target = _field_error(req.target_profile, metrics)
    penalties = compute_penalties(sim_on_target, req.target_profile, metrics, req, coils)
    total = base_mse +
        req.penalties.field * penalties.field +
        req.penalties.geometry * penalties.geometry +
        req.penalties.safety * penalties.safety +
        req.penalties.conductor * penalties.conductor
    return total, base_mse, penalties
end

# Exports
export CoilDescriptor, SimulationMetrics, PenaltyBreakdown, objective_function

end # module


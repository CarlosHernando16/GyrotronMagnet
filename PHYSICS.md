# Electromagnetic Formulation for Gyrotron Magnet Design

## Governing Equations

Magnetostatic problem with magnetic vector potential A:
- ∇×(ν∇×A) = J in Ω
- A×n = 0 on ∂Ω (or appropriate BC)
- ν = 1/μ (reluctivity)
- J = current density from coils

Magnetic field: B = ∇×A

## Weak Form

Find A ∈ V such that:
∫_Ω (ν∇×A)·(∇×v) dΩ = ∫_Ω J·v dΩ, ∀v ∈ V

## Optimization Problem

min_θ f(θ) = ||B(θ) - B_target||²_L2
subject to:
- I_min ≤ I_i ≤ I_max (current bounds)
- r_min ≤ r_i ≤ r_max (geometric bounds)
- Manufacturing constraints

where θ = [r₁, z₁, I₁, ..., rₙ, zₙ, Iₙ]

## Physical Parameters

- Target field: B_z along axis (0, 0, z)
- Typical values: 1-10 Tesla
- Coil currents: 100-10000 A
- Geometry: cylindrical symmetry (2D axisymmetric)

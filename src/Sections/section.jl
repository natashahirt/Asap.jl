abstract type AbstractCrossSection end

"""
    Section(A, E, G, Ix, Iy, J, ρ=1.0u"kg/m^3"; Ay=Inf*u"m^2", Az=Inf*u"m^2")
    Section(mat::Material, A, Ix, Iy, J; Ay=Inf*u"m^2", Az=Inf*u"m^2")

A cross section assigned to an element.

# Fields
- `A` Area [L²]
- `E` Modulus of Elasticity [Pressure]
- `G` Shear Modulus [Pressure]  
- `Ix` Strong axis moment of inertia [L⁴]
- `Iy` Weak axis moment of inertia [L⁴]
- `J` Torsional constant [L⁴]
- `ρ` Density [M/L³]
- `Ay` Effective shear area in Y direction [L²] (default: Inf for Bernoulli-Euler)
- `Az` Effective shear area in Z direction [L²] (default: Inf for Bernoulli-Euler)

# Beam Formulation
- **Bernoulli-Euler** (default): Set `Ay=Inf, Az=Inf` - shear-rigid, good for slender beams
- **Timoshenko**: Set finite `Ay, Az` - shear-flexible, better for stocky beams

For rectangular sections: `Ay ≈ Az ≈ 5/6 * A`
For I-sections: `Ay ≈ web_area`, `Az ≈ 5/6 * flange_area`

# Examples
```julia
# Bernoulli-Euler beam (default, shear-rigid)
Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")

# Timoshenko beam (shear-flexible)
Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4"; 
        Ay=0.0083u"m^2", Az=0.0083u"m^2")
```
"""
struct Section <: AbstractCrossSection
    A::Area                 # Area
    E::Pressure             # Young's modulus
    G::Pressure             # Shear modulus
    Ix::MomentOfInertia     # Strong axis I
    Iy::MomentOfInertia     # Weak axis I
    J::MomentOfInertia      # Torsional constant
    ρ::Density              # Density
    Ay::Area                # Effective shear area Y (Inf m² for Bernoulli-Euler)
    Az::Area                # Effective shear area Z (Inf m² for Bernoulli-Euler)

    function Section(
        A::Quantity,
        E::Quantity,
        G::Quantity,
        Ix::Quantity,
        Iy::Quantity,
        J::Quantity,
        ρ::Quantity = 1.0u"kg/m^3";
        Ay::Union{Quantity, Real} = Inf * u"m^2",
        Az::Union{Quantity, Real} = Inf * u"m^2"
    )
        A_si  = uconvert(u"m^2", A)
        E_si  = uconvert(u"Pa", E)
        G_si  = uconvert(u"Pa", G)
        Ix_si = uconvert(u"m^4", Ix)
        Iy_si = uconvert(u"m^4", Iy)
        J_si  = uconvert(u"m^4", J)
        ρ_si  = uconvert(u"kg/m^3", ρ)
        Ay_si = Ay isa Quantity ? uconvert(u"m^2", Ay) : Ay * u"m^2"
        Az_si = Az isa Quantity ? uconvert(u"m^2", Az) : Az * u"m^2"
        return new(A_si, E_si, G_si, Ix_si, Iy_si, J_si, ρ_si, Ay_si, Az_si)
    end

    function Section(mat::Material, A::Quantity, Ix::Quantity, Iy::Quantity, J::Quantity;
                     Ay::Union{Quantity, Real} = Inf * u"m^2",
                     Az::Union{Quantity, Real} = Inf * u"m^2")
        A_si  = uconvert(u"m^2", A)
        Ix_si = uconvert(u"m^4", Ix)
        Iy_si = uconvert(u"m^4", Iy)
        J_si  = uconvert(u"m^4", J)
        Ay_si = Ay isa Quantity ? uconvert(u"m^2", Ay) : Ay * u"m^2"
        Az_si = Az isa Quantity ? uconvert(u"m^2", Az) : Az * u"m^2"
        return new(A_si, mat.E, mat.G, Ix_si, Iy_si, J_si, mat.ρ, Ay_si, Az_si)
    end
end

"""
Check if section uses Timoshenko (shear-flexible) formulation.
"""
is_timoshenko(sec::Section) = isfinite(sec.Ay) && isfinite(sec.Az)

"""
Check if section uses Bernoulli-Euler (shear-rigid) formulation.
"""
is_bernoulli_euler(sec::Section) = !is_timoshenko(sec)

"""
    TrussSection(A, E, ρ=1.0u"kg/m^3")
    TrussSection(mat::Material, A)

A cross section assigned to a truss element.

# Fields
- `A` Area [L²]
- `E` Modulus of Elasticity [Pressure]
- `ρ` Density [M/L³]

# Examples
```julia
TrussSection(0.01u"m^2", 200u"GPa", 7850u"kg/m^3")
TrussSection(mat, 0.01u"m^2")
```
"""
struct TrussSection <: AbstractCrossSection
    A::Area                 # Area
    E::Pressure             # Young's modulus
    ρ::Density              # Density

    function TrussSection(
        A::Quantity,
        E::Quantity,
        ρ::Quantity = 1.0u"kg/m^3"
    )
        # Convert all to base SI units
        A_si = uconvert(u"m^2", A)
        E_si = uconvert(u"Pa", E)
        ρ_si = uconvert(u"kg/m^3", ρ)
        return new(A_si, E_si, ρ_si)
    end

    function TrussSection(mat::Material, A::Quantity)
        A_si = uconvert(u"m^2", A)
        return new(A_si, mat.E, mat.ρ)
    end
end

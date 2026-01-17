abstract type AbstractSection end

"""
    Section(A, E, G, Ix, Iy, J, ρ = 1.0u"kg/m^3")
    Section(mat::Material, A, Ix, Iy, J)

A cross section assigned to an element.

# Fields
- `A` Area [Distance²] - Quantity
- `E` Modulus of Elasticity [Force/Distance²] - Quantity{Pressure}
- `G` Shear Modulus [Force/Distance²] - Quantity{Pressure}
- `Ix` Nominal strong moment of inertia [Distance⁴] - Quantity
- `Iy` Nominal weak moment of inertia [Distance⁴] - Quantity
- `J` Torsional constant [Distance⁴] - Quantity
- `ρ` Density [Mass/Distance³] - Quantity{Density}

# Examples
```julia
Section(0.01u"m^2", 200u"GPa", 77u"GPa", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4", 7850u"kg/m^3")
Section(mat, 0.01u"m^2", 1e-4u"m^4", 1e-4u"m^4", 1e-4u"m^4")
```
"""
struct Section <: AbstractSection
    A::Area # area
    E::QuantityPressure # young's modulus
    G::QuantityPressure # shear modulus
    Ix::MomentOfInertia # strong axis I
    Iy::MomentOfInertia # weak axis I
    J::MomentOfInertia # torsional constant
    ρ::QuantityDensity # density

    function Section(
        A::Quantity,
        E::Quantity,
        G::Quantity,
        Ix::Quantity,
        Iy::Quantity,
        J::Quantity,
        ρ::Quantity = 1.0u"kg/m^3"
    )
        # Convert all to base SI units
        A_si = uconvert(u"m^2", A)
        E_si = uconvert(u"Pa", E)
        G_si = uconvert(u"Pa", G)
        Ix_si = uconvert(u"m^4", Ix)
        Iy_si = uconvert(u"m^4", Iy)
        J_si = uconvert(u"m^4", J)
        ρ_si = uconvert(u"kg/m^3", ρ)
        return new(A_si, E_si, G_si, Ix_si, Iy_si, J_si, ρ_si)
    end

    function Section(mat::Material, A::Quantity, Ix::Quantity, Iy::Quantity, J::Quantity)
        A_si = uconvert(u"m^2", A)
        Ix_si = uconvert(u"m^4", Ix)
        Iy_si = uconvert(u"m^4", Iy)
        J_si = uconvert(u"m^4", J)
        return new(A_si, mat.E, mat.G, Ix_si, Iy_si, J_si, mat.ρ)
    end
end

"""
    TrussSection(A, E, ρ = 1.0u"kg/m^3")
    TrussSection(mat::Material, A)

A cross section assigned to a truss element.

# Fields
- `A` Area [Distance²] - Quantity
- `E` Modulus of Elasticity [Force/Distance²] - Quantity{Pressure}
- `ρ` Density [Mass/Distance³] - Quantity{Density}

# Examples
```julia
TrussSection(0.01u"m^2", 200u"GPa", 7850u"kg/m^3")
TrussSection(mat, 0.01u"m^2")
```
"""
struct TrussSection <: AbstractSection
    A::Area # area
    E::QuantityPressure # young's modulus
    ρ::QuantityDensity # density

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

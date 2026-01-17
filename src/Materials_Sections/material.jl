"""
    Material(E, G, ρ, ν)

Define a structural material.

# Fields
- `E` Modulus of Elasticity [Force/distance²] - Quantity{Pressure}
- `G` Shear Modulus [Force/distance²] - Quantity{Pressure}
- `ρ` Density [Mass/distance³] - Quantity{Density}
- `ν` Poisson's Ratio [unitless] - Float64

# Examples
```julia
Material(200u"GPa", 77u"GPa", 7850u"kg/m^3", 0.3)
```
"""
struct Material
    E::QuantityPressure #young's modulus
    G::QuantityPressure #shear modulus
    ρ::QuantityDensity #density
    ν::Float64 #poisson's ratio

    function Material(E::Quantity, G::Quantity, ρ::Quantity, ν::Real)
        # Convert all to base SI units
        E_si = uconvert(u"Pa", E)
        G_si = uconvert(u"Pa", G)
        ρ_si = uconvert(u"kg/m^3", ρ)
        return new(E_si, G_si, ρ_si, Float64(ν))
    end
end

# Constants with explicit units
const Steel_Nmm = Material(200e3u"N/mm^2", 77e3u"N/mm^2", 7850.0u"kg/m^3", 0.3)
const Steel_kNm = Material(200e6u"N/m^2", 77e6u"N/m^2", 7850.0u"kg/m^3", 0.3)
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
    E::Pressure   # Young's modulus
    G::Pressure   # Shear modulus
    ρ::Density    # Density
    ν::Float64 #poisson's ratio

    function Material(E::Quantity, G::Quantity, ρ::Quantity, ν::Real)
        # Convert all to base SI units
        E_si = uconvert(u"Pa", E)
        G_si = uconvert(u"Pa", G)
        ρ_si = uconvert(u"kg/m^3", ρ)
        return new(E_si, G_si, ρ_si, Float64(ν))
    end
end

# Note: Material presets removed. Use StructuralSizer materials with to_asap_section().
# Example: to_asap_section(section, A992_Steel) converts to Asap.Section

# =============================================================================
# ShellMaterial - Simplified material for shell elements
# =============================================================================

"""
    ShellMaterial(; E, ν, ρ=0.0, name=:material)

Define a material for shell elements. Simpler than `Material` since shells
derive shear modulus from E and ν.

# Arguments (keyword)
- `E`: Young's modulus (e.g., `30u"GPa"`, `4000u"ksi"`)
- `ν`: Poisson's ratio (0 < ν < 0.5)
- `ρ`: Mass density (default: `0.0u"kg/m^3"`)
- `name`: Optional identifier symbol (default: `:material`)

# Examples
```julia
# Concrete
concrete = ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)

# Steel deck
steel = ShellMaterial(E=200u"GPa", ν=0.3, ρ=7850u"kg/m^3", name=:steel)

# Iteration-friendly usage
for mat in [concrete, steel]
    shells = [ShellTri3(nodes, 0.2u"m", mat) for nodes in mesh]
    # ...
end
```
"""
struct ShellMaterial
    E::Float64      # Young's modulus [Pa]
    ν::Float64      # Poisson's ratio [-]
    ρ::Float64      # Density [kg/m³]
    name::Symbol    # Identifier

    function ShellMaterial(; 
        E::Quantity, 
        ν::Real, 
        ρ::Quantity = 0.0u"kg/m^3",
        name::Symbol = :material
    )
        @assert 0.0 < ν < 0.5 "Poisson's ratio must be in (0, 0.5)"
        
        E_si = ustrip(u"Pa", E)
        ρ_si = ustrip(u"kg/m^3", ρ)
        
        new(E_si, Float64(ν), ρ_si, name)
    end
end

# Note: ShellMaterial presets removed. Create materials inline or use StructuralSizer materials.
# Example: ShellMaterial(E=30u"GPa", ν=0.2, ρ=2400u"kg/m^3", name=:concrete)
# Unit type definitions and aliases for Asap.
# Provides type-safe unit handling using Unitful.jl.

using Unitful
using Unitful: lbf  # Built-in Unitful unit
using StructuralBase.StructuralUnits: kip, ksi  # Custom units from StructuralBase

# =============================================================================
# Base Unit Types (concrete type unions for input validation)
# =============================================================================

"""Distance units: m, mm, cm, ft, inch."""
const Distance = Union{
    typeof(1.0u"m"), typeof(1.0u"mm"), typeof(1.0u"cm"),
    typeof(1.0u"ft"), typeof(1.0u"inch")
}

"""Force units: N, kN, MN, lbf, kip."""
const Force = Union{
    typeof(1.0u"N"), typeof(1.0u"kN"), typeof(1.0u"MN"),
    typeof(1.0lbf), typeof(1.0kip)
}

"""Pressure/stress units: Pa, kPa, MPa, GPa, psi, ksi."""
const Pressure = Union{
    typeof(1.0u"Pa"), typeof(1.0u"kPa"), typeof(1.0u"MPa"),
    typeof(1.0u"GPa"), typeof(1.0u"psi"), typeof(1.0ksi)
}

"""Density units: kg/m³, g/cm³, lb/ft³."""
const Density = Union{
    typeof(1.0u"kg/m^3"), typeof(1.0u"g/cm^3"), typeof(1.0u"lb/ft^3")
}

"""Acceleration units: m/s², ft/s²."""
const Acceleration = Union{typeof(1.0u"m/s^2"), typeof(1.0u"ft/s^2")}

# =============================================================================
# Derived Unit Types (SI storage types)
# =============================================================================

"""Length in m (internal storage)."""
const Length = typeof(1.0u"m")

"""Area in m² (internal storage)."""
const Area = typeof(1.0u"m^2")

"""Volume in m³ (internal storage)."""
const Volume = typeof(1.0u"m^3")

"""Moment of inertia in m⁴ (internal storage)."""
const MomentOfInertia = typeof(1.0u"m^4")

"""Distributed load in N/m (internal storage)."""
const ForcePerLength = typeof(1.0u"N/m")

"""Moment in N*m (internal storage)."""
const Moment = typeof(1.0u"N*m")

# =============================================================================
# Quantity Type Aliases (dimension-based, any compatible unit)
# =============================================================================

"""Any quantity with distance dimension (𝐋)."""
const QuantityDistance = Quantity{T, Unitful.𝐋, U} where {T<:Real, U}

"""Any quantity with force dimension (𝐌𝐋𝐓⁻²)."""
const QuantityForce = Quantity{T, Unitful.𝐌*Unitful.𝐋*Unitful.𝐓^-2, U} where {T<:Real, U}

"""Any quantity with pressure dimension (𝐌𝐋⁻¹𝐓⁻²)."""
const QuantityPressure = Quantity{T, Unitful.𝐌*Unitful.𝐋^-1*Unitful.𝐓^-2, U} where {T<:Real, U}

"""Any quantity with density dimension (𝐌𝐋⁻³)."""
const QuantityDensity = Quantity{T, Unitful.𝐌*Unitful.𝐋^-3, U} where {T<:Real, U}

"""Any quantity with acceleration dimension (𝐋𝐓⁻²)."""
const QuantityAcceleration = Quantity{T, Unitful.𝐋*Unitful.𝐓^-2, U} where {T<:Real, U}

# =============================================================================
# Physical Constants
# =============================================================================

"""Standard gravity acceleration."""
const STANDARD_GRAVITY = 9.80665u"m/s^2"

"""Standard gravity as Float64 (for backward compatibility)."""
const STANDARD_GRAVITY_F64 = 9.80665

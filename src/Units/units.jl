# =============================================================================
# Asap Unit System
# =============================================================================
# Canonical source of structural engineering units for the ecosystem.
# All other packages (StructuralSizer, StructuralSynthesizer) import from here.
#
# Contents:
# - US Customary units (kip, ksi, psf, ksf, pcf)
# - Dimension-based type aliases (Length, Area, Pressure, Force, Moment, etc.)
# - Concrete type aliases for struct fields (LengthQuantity, AreaQuantity, etc.)
# - Unit conversion helpers (to_ksi, to_inches, to_meters, etc.)
# - Catalog parsing utilities (asfloat, maybe_asfloat)

using Unitful
using Unitful: lbf  # Built-in Unitful unit

# =============================================================================
# US Customary Units
# =============================================================================

"""Kilopound-force: 1 kip = 1000 lbf = 4448.22 N"""
Unitful.@unit kip "kip" Kip 4448.2216152605u"N" false

"""Kips per square inch: 1 ksi = 1000 psi = 6.895 MPa"""
Unitful.@unit ksi "ksi" KipPerSquareInch 6.894757e6u"Pa" false

"""Pounds per square foot: 1 psf = 1 lbf/ft² = 47.88 Pa"""
Unitful.@unit psf "psf" PoundPerSquareFoot 47.88025898u"Pa" false

"""Kips per square foot: 1 ksf = 1000 psf"""
Unitful.@unit ksf "ksf" KipPerSquareFoot 47880.25898u"Pa" false

"""Pounds per cubic foot: 1 pcf = 1 lb/ft³"""
Unitful.@unit pcf "pcf" PoundPerCubicFoot 16.01846337u"kg/m^3" false

# =============================================================================
# Physical Constants
# =============================================================================

"""Standard gravity acceleration."""
const GRAVITY = 9.80665u"m/s^2"

"""Standard gravity as Float64 (for backward compatibility)."""
const STANDARD_GRAVITY = GRAVITY
const STANDARD_GRAVITY_F64 = 9.80665

# =============================================================================
# Dimension-Based Type Aliases
# =============================================================================
# These enable cleaner type annotations: `f(x::Length)` instead of complex Quantity types

"""Length quantity (m, ft, inch, etc.)"""
const Length = Unitful.Quantity{T, Unitful.𝐋, U} where {T<:Real, U}

"""Area quantity (m², ft², inch², etc.)"""
const Area = Unitful.Quantity{T, Unitful.𝐋^2, U} where {T<:Real, U}

"""Volume or section modulus quantity (m³, ft³, inch³, etc.)"""
const Volume = Unitful.Quantity{T, Unitful.𝐋^3, U} where {T<:Real, U}

"""Second moment of area (m⁴, ft⁴, inch⁴, etc.) - used for beam bending I = ∫y²dA"""
const SecondMomentOfArea = Unitful.Quantity{T, Unitful.𝐋^4, U} where {T<:Real, U}

"""Torsional constant (same L⁴ dimension as SecondMomentOfArea, but different concept)"""
const TorsionalConstant = SecondMomentOfArea

"""Alias for SecondMomentOfArea"""
const MomentOfInertia = SecondMomentOfArea

"""Warping constant quantity (inch⁶, etc.)"""
const WarpingConstant = Unitful.Quantity{T, Unitful.𝐋^6, U} where {T<:Real, U}

"""Pressure/stress quantity (Pa, psi, ksi, etc.)"""
const Pressure = Unitful.Quantity{T, Unitful.𝐌*Unitful.𝐋^-1*Unitful.𝐓^-2, U} where {T<:Real, U}

"""Force quantity (N, lbf, kip, etc.)"""
const Force = Unitful.Quantity{T, Unitful.𝐌*Unitful.𝐋*Unitful.𝐓^-2, U} where {T<:Real, U}

"""Moment/torque quantity (N·m, kip·ft, lb·in, etc.) - same dimension as Energy"""
const Moment = Unitful.Quantity{T, Unitful.𝐌*Unitful.𝐋^2*Unitful.𝐓^-2, U} where {T<:Real, U}

"""Torque quantity - alias for Moment (same physical dimension: M·L²·T⁻²)"""
const Torque = Moment

"""Linear load quantity (N/m, kip/ft, plf, etc.) - Force per unit length"""
const LinearLoad = Unitful.Quantity{T, Unitful.𝐌*Unitful.𝐓^-2, U} where {T<:Real, U}

"""Area load quantity (Pa, psf, ksf, etc.) - Force per unit area (same as Pressure)"""
const AreaLoadQuantity = Pressure

"""Density quantity (kg/m³, pcf, etc.)"""
const Density = Unitful.Quantity{T, Unitful.𝐌*Unitful.𝐋^-3, U} where {T<:Real, U}

"""Acceleration quantity (m/s², ft/s², etc.)"""
const Acceleration = Unitful.Quantity{T, Unitful.𝐋*Unitful.𝐓^-2, U} where {T<:Real, U}

# Legacy aliases for Asap internal use
const QuantityDistance = Length
const QuantityForce = Force
const QuantityPressure = Pressure
const QuantityDensity = Density
const QuantityAcceleration = Acceleration

# =============================================================================
# Concrete Type Aliases (for struct fields with Float64 precision)
# =============================================================================
# Use these when you need concrete types for struct fields.
# The abstract types above (Length, Area, etc.) are for function signatures.

"""Concrete Length type in meters (Float64 precision)."""
const LengthQuantity = typeof(1.0u"m")

"""Concrete Area type in square meters (Float64 precision)."""
const AreaQuantity = typeof(1.0u"m^2")

"""Concrete Volume type in cubic meters (Float64 precision)."""
const VolumeQuantity = typeof(1.0u"m^3")

"""Concrete Pressure type in Pascals (Float64 precision)."""
const PressureQuantity = typeof(1.0u"Pa")

"""Concrete Force type in Newtons (Float64 precision)."""
const ForceQuantity = typeof(1.0u"N")

"""Concrete Moment type in Newton-meters (Float64 precision)."""
const MomentQuantity = typeof(1.0u"N*m")

"""Distributed load in N/m (internal storage)."""
const ForcePerLength = typeof(1.0u"N/m")

# =============================================================================
# Concrete Type Unions (for input validation)
# =============================================================================

"""Distance units: m, mm, cm, ft, inch."""
const Distance = Union{
    typeof(1.0u"m"), typeof(1.0u"mm"), typeof(1.0u"cm"),
    typeof(1.0u"ft"), typeof(1.0u"inch")
}

# =============================================================================
# Unit Conversion Helpers - US Customary
# =============================================================================

"""
    to_inches(x) -> Float64

Convert a length to inches (stripped of units).
If already a Real number, assumes it's in inches and returns as Float64.
"""
to_inches(x::Length) = Float64(ustrip(u"inch", x))
to_inches(x::Real) = Float64(x)

"""
    to_sqinches(x) -> Float64

Convert an area to square inches (stripped of units).
If already a Real number, assumes it's in in² and returns as Float64.
"""
to_sqinches(x::Area) = Float64(ustrip(u"inch^2", x))
to_sqinches(x::Real) = Float64(x)

"""
    to_ksi(x) -> Float64

Convert a pressure to ksi (stripped of units).
If already a Real number, assumes it's in ksi and returns as Float64.
"""
to_ksi(x::Pressure) = Float64(ustrip(ksi, x))
to_ksi(x::Real) = Float64(x)

"""
    to_kip(x) -> Float64

Convert a force to kip (stripped of units).
If already a Real number, assumes it's in kip and returns as Float64.
"""
to_kip(x::Force) = Float64(ustrip(kip, x))
to_kip(x::Real) = Float64(x)

"""
    to_kipft(x) -> Float64

Convert a moment to kip-ft (stripped of units).
If already a Real number, assumes it's in kip-ft and returns as Float64.
"""
to_kipft(x::Moment) = Float64(ustrip(kip*u"ft", x))
to_kipft(x::Real) = Float64(x)

# =============================================================================
# Unit Conversion Helpers - SI
# =============================================================================

"""
    to_meters(x) -> Float64

Convert a length to meters (stripped of units).
If already a Real number, assumes it's in meters and returns as Float64.
"""
to_meters(x::Length) = Float64(ustrip(u"m", x))
to_meters(x::Quantity) = ustrip(u"m", x)  # Fallback for any Quantity
to_meters(x::Real) = Float64(x)

"""Convert to m² and strip units."""
to_meters_squared(x::Quantity) = ustrip(u"m^2", x)
to_meters_squared(x::Real) = Float64(x)

"""Convert to m⁴ and strip units."""
to_meters_fourth(x::Quantity) = ustrip(u"m^4", x)
to_meters_fourth(x::Real) = Float64(x)

"""
    to_pascals(x) -> Float64

Convert a pressure to Pascals (stripped of units).
If already a Real number, assumes it's in Pa and returns as Float64.
"""
to_pascals(x::Pressure) = Float64(ustrip(u"Pa", x))
to_pascals(x::Real) = Float64(x)

"""
    to_newtons(x) -> Float64

Convert a force to Newtons (stripped of units).
If already a Real number, assumes it's in N and returns as Float64.
"""
to_newtons(x::Force) = Float64(ustrip(u"N", x))
to_newtons(x::Quantity) = ustrip(u"N", x)  # Fallback
to_newtons(x::Real) = Float64(x)

"""
    to_newton_meters(x) -> Float64

Convert a moment to Newton-meters (stripped of units).
If already a Real number, assumes it's in N·m and returns as Float64.
"""
to_newton_meters(x::Moment) = Float64(ustrip(u"N*m", x))
to_newton_meters(x::Quantity) = ustrip(u"N*m", x)  # Fallback
to_newton_meters(x::Real) = Float64(x)

"""
    to_newtons_per_meter(x) -> Float64

Convert a linear load to N/m (stripped of units).
If already a Real number, assumes it's in N/m and returns as Float64.
"""
to_newtons_per_meter(x::LinearLoad) = Float64(ustrip(u"N/m", x))
to_newtons_per_meter(x::Quantity) = ustrip(u"N/m", x)  # Fallback
to_newtons_per_meter(x::Real) = Float64(x)

"""Convert to kg/m³ and strip units."""
to_kg_per_m3(x::Quantity) = ustrip(u"kg/m^3", x)
to_kg_per_m3(x::Real) = Float64(x)

"""Convert to m/s² and strip units."""
to_m_per_s2(x::Quantity) = ustrip(u"m/s^2", x)
to_m_per_s2(x::Real) = Float64(x)

# =============================================================================
# Displacement/Reaction Vector Converters
# =============================================================================

"""
Convert a displacement vector to Float64.
DOFs 1-3 are translations (m), DOFs 4-6 are rotations (rad).
"""
function to_displacement_vec(v::Vector{Quantity})
    result = Vector{Float64}(undef, length(v))
    for (i, d) in enumerate(v)
        result[i] = i <= 3 ? ustrip(u"m", d) : ustrip(d)
    end
    return result
end

"""
Convert a reaction vector to Float64.
DOFs 1-3 are forces (N), DOFs 4-6 are moments (N*m).
"""
function to_reaction_vec(v::Vector{Quantity})
    result = Vector{Float64}(undef, length(v))
    for (i, r) in enumerate(v)
        result[i] = i <= 3 ? ustrip(u"N", r) : ustrip(u"N*m", r)
    end
    return result
end

# =============================================================================
# CSV/Catalog Parsing Utilities
# =============================================================================

"""Convert a value to Float64 (for catalog parsing)."""
asfloat(x::Real) = Float64(x)
asfloat(x::AbstractString) = parse(Float64, x)
asfloat(x) = throw(ArgumentError("Cannot parse numeric value from $(typeof(x)) = $(repr(x))"))

"""Convert a value to Float64 or nothing if missing/invalid (for optional catalog fields)."""
function maybe_asfloat(x)
    ismissing(x) && return nothing
    if x isa AbstractString
        sx = strip(x)
        (sx == "–" || sx == "-" || sx == "—" || isempty(sx)) && return nothing
    end
    return asfloat(x)
end

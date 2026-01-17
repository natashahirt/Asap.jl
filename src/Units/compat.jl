using Unitful

"""
Compatibility layer for Unitful migration.

Provides helper functions for converting between Unitful and Float64,
and for handling backward compatibility during the migration.
"""

# =============================================================================
# Unit Conversion Helpers
# =============================================================================
# Note: ustrip(unit, x) automatically converts x to the target unit before
# stripping, so explicit uconvert calls are unnecessary.

"""Convert a distance quantity to meters and strip units."""
to_meters(x::Quantity) = ustrip(u"m", x)
to_meters(x::Real) = Float64(x)

"""Convert an area quantity to m² and strip units."""
to_meters_squared(x::Quantity) = ustrip(u"m^2", x)
to_meters_squared(x::Real) = Float64(x)

"""Convert a moment of inertia quantity to m⁴ and strip units."""
to_meters_fourth(x::Quantity) = ustrip(u"m^4", x)
to_meters_fourth(x::Real) = Float64(x)

"""Convert a force quantity to Newtons and strip units."""
to_newtons(x::Quantity) = ustrip(u"N", x)
to_newtons(x::Real) = Float64(x)

"""Convert a pressure/stress quantity to Pascals and strip units."""
to_pascals(x::Quantity) = ustrip(u"Pa", x)
to_pascals(x::Real) = Float64(x)

"""Convert a density quantity to kg/m³ and strip units."""
to_kg_per_m3(x::Quantity) = ustrip(u"kg/m^3", x)
to_kg_per_m3(x::Real) = Float64(x)

"""Convert an acceleration quantity to m/s² and strip units."""
to_m_per_s2(x::Quantity) = ustrip(u"m/s^2", x)
to_m_per_s2(x::Real) = Float64(x)

# =============================================================================
# Vector Conversion Helpers
# =============================================================================

"""Convert a vector of distance quantities to meters."""
to_meters_vec(v::Vector{<:Quantity}) = [to_meters(x) for x in v]
to_meters_vec(v::Vector{<:Real}) = Float64.(v)

"""Convert a vector of force quantities to Newtons."""
to_newtons_vec(v::Vector{<:Quantity}) = [to_newtons(x) for x in v]
to_newtons_vec(v::Vector{<:Real}) = Float64.(v)

"""Convert a vector of pressure quantities to Pascals."""
to_pascals_vec(v::Vector{<:Quantity}) = [to_pascals(x) for x in v]
to_pascals_vec(v::Vector{<:Real}) = Float64.(v)

# =============================================================================
# Promotion Helpers
# =============================================================================

"""Promote a Real value to a Quantity with the given unit."""
promote_to_quantity(x::Real, unit::Unitful.Units) = x * unit
promote_to_quantity(x::Quantity, unit::Unitful.Units) = uconvert(unit, x)

"""Promote a Real to distance (assumes meters)."""
promote_to_distance(x::Real) = x * u"m"
promote_to_distance(x::Quantity) = uconvert(u"m", x)

"""Promote a Real to force (assumes Newtons)."""
promote_to_force(x::Real) = x * u"N"
promote_to_force(x::Quantity) = uconvert(u"N", x)

"""Promote a Real to pressure (assumes Pascals)."""
promote_to_pressure(x::Real) = x * u"Pa"
promote_to_pressure(x::Quantity) = uconvert(u"Pa", x)

"""Promote a Real to density (assumes kg/m³)."""
promote_to_density(x::Real) = x * u"kg/m^3"
promote_to_density(x::Quantity) = uconvert(u"kg/m^3", x)

"""Promote a Real to acceleration (assumes m/s²)."""
promote_to_acceleration(x::Real) = x * u"m/s^2"
promote_to_acceleration(x::Quantity) = uconvert(u"m/s^2", x)

"""Convert a moment quantity to N*m and strip units."""
to_newton_meters(x::Quantity) = ustrip(u"N*m", x)
to_newton_meters(x::Real) = Float64(x)

"""Convert a force per length quantity to N/m and strip units."""
to_newtons_per_meter(x::Quantity) = ustrip(u"N/m", x)
to_newtons_per_meter(x::Real) = Float64(x)

# =============================================================================
# Vector Promotion Helpers
# =============================================================================

"""Promote a vector to distance quantities (assumes meters)."""
promote_to_distance_vec(v::Vector{<:Real}) = [promote_to_distance(x) for x in v]
promote_to_distance_vec(v::Vector{<:Quantity}) = [uconvert(u"m", x) for x in v]

"""Promote a vector to force quantities (assumes Newtons)."""
promote_to_force_vec(v::Vector{<:Real}) = [promote_to_force(x) for x in v]
promote_to_force_vec(v::Vector{<:Quantity}) = [uconvert(u"N", x) for x in v]

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

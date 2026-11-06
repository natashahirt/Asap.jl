
# =============================================================================
# Helper: rotate a 3-vector from global to local with beam sign convention
# =============================================================================

"""
    _rotate_to_local(Rv, g1, g2, g3)

Rotate global vector (g1,g2,g3) into the local frame using the 3×3 rotation
sub-block `Rv = @view element.R[1:3,1:3]`.  Returns `(wx, wy, wz)` with the
beam sign convention (y,z flipped).
"""
@inline function _rotate_to_local(Rv, g1, g2, g3)
    @inbounds begin
        wx =  Rv[1,1]*g1 + Rv[1,2]*g2 + Rv[1,3]*g3
        wy = -(Rv[2,1]*g1 + Rv[2,2]*g2 + Rv[2,3]*g3)
        wz = -(Rv[3,1]*g1 + Rv[3,2]*g2 + Rv[3,3]*g3)
    end
    return wx, wy, wz
end

# =============================================================================
# ElementProps — pre-stripped element properties for hot-path loops
# =============================================================================

"""
    ElementProps

Pre-stripped (unitless Float64) element properties, computed once per element
and passed into `accumulate_force!` / `accumulatedisp!` to avoid repeated
`ustrip` calls inside inner loops.
"""
struct ElementProps
    Rv::SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false}
    L::Float64
    E::Float64
    Ix::Float64   # strong-axis I (I₂ / Iy in Asap convention)
    Iy::Float64   # weak-axis I  (I₃ / Iz in Asap convention)
end

"""Pre-strip element properties once for hot-path force/displacement loops."""
@inline function ElementProps(element::AbstractElement)
    ElementProps(
        (@view element.R[1:3, 1:3]),
        ustrip(u"m", element.length),
        ustrip(u"Pa", element.section.E),
        ustrip(u"m^4", element.section.Ix),
        ustrip(u"m^4", element.section.Iy),
    )
end

# =============================================================================
# _fill_uglobal! — fill 12-vec with nodal displacements (avoids to_displacement_vec + concat)
# =============================================================================

"""
    _fill_uglobal!(buf, element)

Fill a 12-element buffer with nodal displacements directly from element node
data.  DOFs 1-3 are translations (m), DOFs 4-6 are rotations (rad).
Avoids two `to_displacement_vec` allocations + `[start; end]` concatenation.
"""
@inline function _fill_uglobal!(buf::Vector{Float64}, element::AbstractElement)
    ds = element.nodeStart.displacement
    de = element.nodeEnd.displacement
    @inbounds begin
        buf[1]  = ustrip(u"m", ds[1])
        buf[2]  = ustrip(u"m", ds[2])
        buf[3]  = ustrip(u"m", ds[3])
        buf[4]  = ustrip(ds[4])
        buf[5]  = ustrip(ds[5])
        buf[6]  = ustrip(ds[6])
        buf[7]  = ustrip(u"m", de[1])
        buf[8]  = ustrip(u"m", de[2])
        buf[9]  = ustrip(u"m", de[3])
        buf[10] = ustrip(de[4])
        buf[11] = ustrip(de[5])
        buf[12] = ustrip(de[6])
    end
    return buf
end


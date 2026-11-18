const etype2DOF = Dict(
    Element{FixedFixed} => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    Element{FreeFree} => [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0],
    Element{FixedFree} => [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    Element{FreeFixed} => [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    Element{Joist} => [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0]
)

"""
    release_end_rotations!(ul, element, L) -> ul

Replace the bending rotations of `element`'s released ends in the local
12-vector `ul` (after the `etype2DOF` mask) with the member's own end slopes.

A released end carries no moment, so the member rotates independently of the
node there. With the slope at a released end set to the value that gives zero
curvature, the Hermite cubic through (v₁, θ₁, v₂, θ₂) is the member's deflected
shape from end motion:

- start released: θ₁ = 3(v₂ - v₁)/(2L) - θ₂/2
- end released:   θ₂ = 3(v₂ - v₁)/(2L) - θ₁/2
- both released:  θ₁ = θ₂ = (v₂ - v₁)/L (straight chord)

In the local xz plane the slope dw/dx is `-ul[5]` (`-ul[11]` at the end).
"""
function release_end_rotations!(ul::AbstractVector{Float64}, element::Element{R}, L::Real) where {R}
    R === FixedFixed && return ul
    start_rel = R === FreeFixed || R === FreeFree || R === Joist
    end_rel = R === FixedFree || R === FreeFree || R === Joist
    # (v₁, v₂, slope index at start, slope index at end, slope sign)
    for (i1, i2, s1, s2, sgn) in ((2, 8, 6, 12, 1.0), (3, 9, 5, 11, -1.0))
        chord = (ul[i2] - ul[i1]) / L
        if start_rel && end_rel
            ul[s1] = sgn * chord
            ul[s2] = sgn * chord
        elseif start_rel
            ul[s1] = sgn * (1.5chord - sgn * ul[s2] / 2)
        elseif end_rel
            ul[s2] = sgn * (1.5chord - sgn * ul[s1] / 2)
        end
    end
    return ul
end

const planarDOFs = Dict(:X => [1, 7],
    :XY => [2, 6, 8, 12],
    :XZ => [3, 5, 9, 11])

"""
Collect all elements that are part of the same continuous member (ie collect shattered elements)
"""
function groupbyid(elements::Vector{<:AbstractElement})
    ichecked = Vector{Int64}()
    indices = Vector{Vector{Int64}}()

    elementids = getproperty.(elements, :elementID)

    for id in elementids
        
        in(id, ichecked) && continue

        igroup = findall(elementids .== id)

        push!(indices, igroup)
        push!(ichecked, id)
    end

    return indices
end
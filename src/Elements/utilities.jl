"""
Extract all elements with a given ID
elements = Vector{Element}()
elements[:outerEdges]
"""
function Base.getindex(elements::Vector{<:AbstractElement}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{<:AbstractElement}, i::Symbol)
    return findall([element.id == i for element in elements])
end

function Base.getindex(elements::Vector{TrussElement}, i::Symbol)
    return [element for element in elements if element.id == i]
end

function Base.findall(elements::Vector{TrussElement}, i::Symbol)
    return findall([element.id == i for element in elements])
end

Base.length(element::T) where {T<:AbstractElement} = norm([to_meters(p) for p in (element.nodeEnd.position .- element.nodeStart.position)]) * u"m"

function length!(element::T) where {T<:AbstractElement}
    L = length(element)
    to_meters(L) > 0 || error("Zero-length element detected between nodes at " *
        "$(element.nodeStart.position) — check for coincident nodes (id=:$(element.id))")
    element.length = L
end


"""
    local_x(element::AbstractElement; unit = true)

Get the local x vector of an element: element.nodeEnd.position - element.nodeStart.position.
`unit = true` gives the normalized vector.
"""
function local_x(element::T; unit = true) where {T<:AbstractElement}
    x_quantities = element.nodeEnd.position .- element.nodeStart.position
    x = [to_meters(p) for p in x_quantities]
    unit ? normalize(x) : x
end

"""
    lcs(element::AbstractElement, Ψ::Float64; tol = 1e-6)

Get the local coordinate system unit vectors of a given element and pitch angle Ψ: [local_x, local_y, local_z]
"""
function lcs(element::T, Ψ::Float64; tol = 1e-6) where {T<:AbstractElement}

    # local x vector
    xvec = local_x(element)
    
    if norm(cross(xvec, globalY)) < tol
        CYx = xvec[2] #cosine to global Y axis
        xvec = CYx * globalY
        yvec = -CYx * globalX * cos(Ψ) + sin(Ψ) * globalZ
        zvec = CYx * globalX * sin(Ψ) + cos(Ψ) * globalZ 
    else
        zbar = normalize(cross(xvec, [0, 1, 0]))
        ybar = normalize(cross(zbar, xvec))

        yvec = cos(Ψ) * ybar + sin(Ψ) * zbar
        zvec = -sin(Ψ) * ybar + cos(Ψ) * zbar
    end

    return [xvec, yvec, zvec]
end

"""
    lcs(element::AbstractElement, Ψ::Float64; tol = 1e-6)

Populate local coordinate system unit vectors of a given element and pitch angle Ψ: [local_x, local_y, local_z]
"""
function lcs!(element::T, Ψ::Float64; tol = 1e-6) where {T<:AbstractElement}

    # local x vector
    xvec = local_x(element)
    
    if norm(cross(xvec, globalY)) < tol
        CYx = xvec[2] #cosine to global Y axis
        xvec = CYx * globalY
        yvec = -CYx * globalX * cos(Ψ) + sin(Ψ) * globalZ
        zvec = CYx * globalX * sin(Ψ) + cos(Ψ) * globalZ 
    else
        zbar = normalize(cross(xvec, [0, 1, 0]))
        ybar = normalize(cross(zbar, xvec))

        yvec = cos(Ψ) * ybar + sin(Ψ) * zbar
        zvec = -sin(Ψ) * ybar + cos(Ψ) * zbar
    end

    element.LCS = [xvec, yvec, zvec]
end


"""
    lcs(xin::Vector{Float64}, Ψ::Float64; tol = 1e-4)

Get the local coordinate system unit vectors of a given local x axis and pitch angle Ψ: [local_x, local_y, local_z]
"""
function lcs(xin::Vector{Float64}, Ψ::Float64; tol = 1e-6)

    xvec = normalize(copy(xin))
    # local x vector
    if norm(cross(xvec, globalY)) < tol
        CYx = xvec[2] #cosine to global Y axis
        xvec = CYx * globalY
        yvec = -CYx * globalX * cos(Ψ) + sin(Ψ) * globalZ
        zvec = CYx * globalX * sin(Ψ) + cos(Ψ) * globalZ 
    else
        zbar = normalize(cross(xvec, [0, 1, 0]))
        ybar = normalize(cross(zbar, xvec))

        yvec = cos(Ψ) * ybar + sin(Ψ) * zbar
        zvec = -sin(Ψ) * ybar + cos(Ψ) * zbar
    end

    return [xvec, yvec, zvec]
end

# A unit local y this close to horizontal has no up/down sense. A truly
# vertical member lands here; a member with a real tilt does not.
const _UP_ECCENTRICITY_HORIZONTAL_TOL = 1e-8

"""
    orient_up_eccentricity(ecc, nodeStart, nodeEnd; roll = π/2) -> NTuple{4, Float64}

Sign a physical upward offset for the eccentricity tuple [`Element`](@ref) stores.

`ecc` is `(e_start_f1, e_end_f1, e_f2, e_f3)`, in metres or length quantities,
the argument order of [`eccentricity_transformation!`](@ref). That transform
is the rigid link `u_slave = u_master + θ × e` with `e = (e_f1, e_f2, e_f3)`
in the local frame, so slot 3 (`e_f2`) is the local-y offset and slot 4
(`e_f3`) is local z.

Pass the physical offset. A positive slot 3 means global up (global +Z).
Do not pass a tuple already signed for this member's local y. The returned
tuple is `ecc` in metres with slot 3 multiplied by `sign(local_y · global Z)`.
Local axes are [`lcs`](@ref) of the chord from `nodeStart` to `nodeEnd` at
`roll`, the same frame `process!` stores on the element. `roll` defaults to
`π/2`, the roll the element constructor stores until a caller sets `Ψ`.
This does not require `process!`.

A vertical member, `|local_y · Z| < 1e-8`, has no up/down seat: slot 3 is
left unchanged and the function does not throw. A zero-length chord has no
local y either, and is likewise returned unchanged.
"""
function orient_up_eccentricity(ecc::NTuple{4}, nodeStart, nodeEnd; roll::Real = pi / 2)
    stored = (
        to_meters(ecc[1]),
        to_meters(ecc[2]),
        to_meters(ecc[3]),
        to_meters(ecc[4]),
    )
    chord = _node_position_m(nodeEnd) .- _node_position_m(nodeStart)
    norm(chord) > 1e-12 || return stored
    local_y = lcs(chord, Float64(roll))[2]
    dot_z = local_y[3]
    abs(dot_z) < _UP_ECCENTRICITY_HORIZONTAL_TOL && return stored
    sense = dot_z > 0 ? 1.0 : -1.0
    return (stored[1], stored[2], stored[3] * sense, stored[4])
end

"""Node position as metres, the units [`lcs`](@ref) consumes."""
function _node_position_m(node)
    pos = node.position
    return Float64[to_meters(pos[1]), to_meters(pos[2]), to_meters(pos[3])]
end

"""
    endpoints(element::AbstractElement)

Extract the start and end points as two vectors.
"""
function endpoints(element::T) where {T<:AbstractElement}
    return [[to_meters(p) for p in element.nodeStart.position], [to_meters(p) for p in element.nodeEnd.position]]
end

"""
    midpoint(element::AbstractElement)

Extract the centerpoint of an element as a vector.
"""
function midpoint(element::T) where {T<:AbstractElement}
    mid_quantities = (element.nodeStart.position .+ element.nodeEnd.position) ./ 2
    return [to_meters(p) for p in mid_quantities]
end

"""
    axial_force(element::TrussElement)

Extract the axial force of an element
"""
function axial_force(element::TrussElement)
    element.forces[2]
end

"""
    axial_force(element::Element)
    
Extract the axial force of an element
"""
function axial_force(element::Element)
    element.forces[7]
end

function nodeids(element::T) where {T<:AbstractElement}
    return [element.nodeStart.nodeID, element.nodeEnd.nodeID]
end
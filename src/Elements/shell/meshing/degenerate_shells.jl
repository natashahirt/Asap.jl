# A triangle with no area makes the shell stiffness undefined. Two corners
# within `tol` are one node. Corners spread along a straight edge stay put,
# and that triangle is removed.

"""Distance between two nodes, in metres."""
function _node_separation(a::Node, b::Node)
    dx = ustrip(u"m", a.position[1]) - ustrip(u"m", b.position[1])
    dy = ustrip(u"m", a.position[2]) - ustrip(u"m", b.position[2])
    dz = ustrip(u"m", a.position[3]) - ustrip(u"m", b.position[3])
    return hypot(dx, dy, dz)
end

"""True when `el` has no area or fewer than three distinct corners."""
function _shell_is_flat(el, min_area::Float64)
    el isa ShellTri3 || return false
    n1, n2, n3 = el.nodes
    (n1 === n2 || n2 === n3 || n3 === n1) && return true
    return el.area < min_area
end

function _degenerate_find!(parent::Dict{UInt64, UInt64}, id::UInt64)
    parent[id] == id && return id
    parent[id] = _degenerate_find!(parent, parent[id])
    return parent[id]
end

"""
Union of corners that sit on each other, taken only from flat triangles.

The representative of a pair is the node that appears first. Nodes that are
further apart than `tol` are not joined, even on the same flat triangle.
"""
function _degenerate_parents(shells, tol::Float64, min_area::Float64)
    parent = Dict{UInt64, UInt64}()
    canon = Dict{UInt64, Node}()
    order = Dict{UInt64, Int}()
    nseen = Ref(0)
    function touch!(node::Node)
        id = objectid(node)
        haskey(parent, id) && return id
        nseen[] += 1
        order[id] = nseen[]
        parent[id] = id
        canon[id] = node
        return id
    end
    for el in shells
        _shell_is_flat(el, min_area) || continue
        ns = el.nodes
        for i in 1:3
            touch!(ns[i])
        end
        for i in 1:2, j in (i + 1):3
            _node_separation(ns[i], ns[j]) <= tol &&
                _degenerate_unite!(parent, order, ns[i], ns[j])
        end
    end
    return parent, canon
end

function _degenerate_unite!(parent, order, a::Node, b::Node)
    ra = _degenerate_find!(parent, objectid(a))
    rb = _degenerate_find!(parent, objectid(b))
    ra == rb && return parent
    if order[ra] > order[rb]
        ra, rb = rb, ra
    end
    parent[rb] = ra
    return parent
end

"""The node `node` was welded to, or `node` when it was not on a flat triangle."""
function _degenerate_rep(node::Node, parent, canon)
    id = objectid(node)
    haskey(parent, id) || return node
    return canon[_degenerate_find!(parent, id)]
end

"""
    collapse_degenerate_shells!(shells; tol=1e-4, min_area=1e-8, alongside=nothing)

Weld nodes that sit on each other, and drop a shell triangle with no area.

Two corners of a flat triangle within `tol` metres are the same point. Every
shell is rewired onto the one that appears first, and a triangle that no
longer has three corners is removed. A triangle that is flat because its
corners lie along a straight edge, metres apart, is removed and those nodes
stay. Welding them would collapse the edge.

`alongside` is a vector parallel to `shells`, such as a panel index. Entries
for removed shells are removed with them.

Call this before the model is assembled. Assembly recomputes element geometry
from the nodes that remain.
"""
function collapse_degenerate_shells!(shells::AbstractVector;
                                     tol::Float64 = 1e-4,
                                     min_area::Float64 = 1e-8,
                                     alongside = nothing)
    alongside === nothing || length(alongside) == length(shells) ||
        throw(ArgumentError("alongside must name every shell"))
    parent, canon = _degenerate_parents(shells, tol, min_area)
    keep = Int[]
    for (i, el) in pairs(shells)
        if el isa ShellTri3 && !isempty(parent)
            rewritten = false
            nodes = ntuple(3) do k
                rep = _degenerate_rep(el.nodes[k], parent, canon)
                rep !== el.nodes[k] && (rewritten = true)
                rep
            end
            if rewritten
                el.nodes = nodes
                el.area = _compute_triangle_area(nodes)
            end
        end
        _shell_is_flat(el, min_area) && continue
        push!(keep, i)
    end
    if length(keep) != length(shells)
        drop = setdiff(eachindex(shells), keep)
        deleteat!(shells, drop)
        alongside === nothing || deleteat!(alongside, drop)
    end
    return shells
end

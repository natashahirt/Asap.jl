"""
Define loads applied to structure/elements
"""
abstract type AbstractLoad end
abstract type NodeLoad <: AbstractLoad end
abstract type ElementLoad{R} <: AbstractLoad end

"""
    NodeForce(node::AbstractNode, value::Vector{Quantity})

A force vector [Fx, Fy, Fz] in the global coordinate system applied to a node.
"""
mutable struct NodeForce <: NodeLoad
    node::Union{Node, TrussNode}
    value::Vector{QuantityForce}
    loadID::Int64
    id::Symbol
    
    function NodeForce(node::AbstractNode, value::Vector{<:Quantity}, id::Symbol = :force)
        @assert length(value) == 3 "load vector must be in R³ (GCS)"
        
        # Convert all to Newtons
        value_si = [uconvert(u"N", v) for v in value]
        force = new(node, value_si, 0, id)
        return force
    end

    function NodeForce(nodes::Vector{<:AbstractNode}, index::Integer, value::Vector{<:Quantity}, id::Symbol = :force)
        return NodeForce(nodes[index], value, id)
    end
end

"""
    NodeMoment(node::Node, value::Vector{Quantity})

A moment vector [Mx, My, Mz] in the global coordinate system applied to a node with rotational DOFs.
"""
mutable struct NodeMoment <: NodeLoad
    node::Node
    value::Vector{Moment}
    loadID::Int64
    id::Symbol
    
    function NodeMoment(node::Node, value::Vector{<:Quantity}, id::Symbol = :moment)
        @assert length(value) == 3 "Moment vector must be in R³ (GCS)"
        
        # Convert all to N*m
        value_si = [uconvert(u"N*m", v) for v in value]
        force = new(node, value_si, 0, id)
        return force
    end
end

"""
    LineLoad(element::Element, value::Vector{Quantity})

A distributed line load [wx, wy, wz] in (force/length) applied along an element in the global coordinate system.
"""
mutable struct LineLoad{R<:Release} <: ElementLoad{R}
    element::FrameElement
    value::Vector{ForcePerLength}
    loadID::Int64
    id::Symbol

    function LineLoad(element::FrameElement{R}, value::Vector{<:Quantity}, id::Symbol = :lineload) where R
        @assert length(value) == 3 "load vector must be in R³ (GCS)"
        
        # Convert all to N/m
        value_si = [uconvert(u"N/m", v) for v in value]
        force = new{R}(element, value_si, 0, id)
        return force
    end
end

"""
    GravityLoad(element::Element, factor::Quantity)

A gravity load (negative global Z) applied along a member. 

Generates distributed load w = element.section.A * element.section.ρ * factor, where factor should be the appropriate acceleration due to gravity.
"""
mutable struct GravityLoad{R<:Release} <: ElementLoad{R}
    element::FrameElement
    factor::QuantityAcceleration
    loadID::Int64
    id::Symbol

    function GravityLoad(element::FrameElement{R}, factor::Quantity, id::Symbol = :gravityload) where R
        # Convert to m/s²
        factor_si = uconvert(u"m/s^2", factor)
        force = new{R}(element, factor_si, 0, id)
        return force
    end
end


"""
    PointLoad(element::Element, position::Float64, value::Vector{Quantity})

A point load [Px, Py, Pz] applied in the global coordinate system at a distance `position` × `element.length` from the starting node.
`position` is dimensionless (fraction 0-1).
"""
mutable struct PointLoad{R<:Release} <: ElementLoad{R}
    element::FrameElement
    position::Float64  # Dimensionless fraction
    value::Vector{QuantityForce}
    loadID::Int64
    id::Symbol

    function PointLoad(element::FrameElement{R}, position::Float64, value::Vector{<:Quantity}, id::Symbol = :pointload) where R
        @assert 0 < position < 1 "position must be ∈ ]0, 1["
        @assert length(value) == 3 "load vector must be in R³ (GCS)"
        
        # Convert all to Newtons
        value_si = [uconvert(u"N", v) for v in value]
        force = new{R}(element, position, value_si, 0, id)
        return force
    end
end
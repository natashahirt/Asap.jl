#=
Element Definitions for Asap.jl
===============================

Core structural element types: frame beams, bridge elements, and truss elements.

Eccentric connection support based on formulations from:
- FinEtoolsFlexStructures.jl by Petr Krysl - MIT License
=#

abstract type AbstractElement end
abstract type FrameElement{R} <: AbstractElement end

# Element release type parameters
abstract type Release end
struct FixedFixed <: Release end
struct FixedFree <: Release end
struct FreeFixed <: Release end
struct FreeFree <: Release end
struct Joist <: Release end

const _ReleaseDict = Dict(
    :fixedfixed => FixedFixed,
    :fixedfree => FixedFree,
    :freefixed => FreeFixed,
    :freefree => FreeFree,
    :joist => Joist
)

"""
    Element(nodeStart, nodeEnd, section, id=:element; release=:fixedfixed, eccentricity=nothing)

Create a frame element with an optional `id` tag.

# Arguments
- `nodeStart::Node`: Start node
- `nodeEnd::Node`: End node  
- `section::Section`: Cross-section properties
- `id::Symbol`: Optional identifier (default: `:element`)

# Keyword Arguments
- `release::Symbol`: End release condition (default: `:fixedfixed`)
- `eccentricity::Union{Nothing, NTuple{4, Quantity}}`: Optional eccentric connection offsets
  - `(e_start_axial, e_end_axial, e_transverse_y, e_transverse_z)`
  - When specified, the beam connects to nodes with these offsets from the node positions
  - Useful for modeling real connections where beam centroids don't align with nodes

# Available Releases
- `:fixedfixed` - all DOFs tied to nodes (default)
- `:fixedfree` - rotational DOFs released at end node
- `:freefixed` - rotational DOFs released at start node  
- `:freefree` - all rotational DOFs released (truss behavior)
- `:joist` - all rotational DOFs except torsion released

# Example
```julia
# Standard beam
Element(node1, node2, section)

# Beam with eccentric connection (offset 0.5m in Y direction)
Element(node1, node2, section; eccentricity=(0.0u"m", 0.0u"m", 0.5u"m", 0.0u"m"))
```
"""
mutable struct Element{R<:Release} <: FrameElement{R}
    section::Section                    # Cross section
    nodeStart::Node                     # Start node
    nodeEnd::Node                       # End node
    elementID::Int64                    # Element index
    globalID::Vector{Int64}             # Element global DOFs
    length::QuantityDistance            # Length of element
    K::Matrix{Float64}                  # Stiffness matrix in GCS
    Q::Vector{Float64}                  # Fixed end forces in GCS
    R::Matrix{Float64}                  # Transformation matrix
    Ψ::Float64                          # Roll angle
    LCS::Vector{Vector{Float64}}        # Local coordinate frame (X, y, z)
    forces::Vector{Float64}             # Elemental forces in LCS
    id::Symbol                          # Optional identifier
    eccentricity::NTuple{4, Float64}    # (e_start_f1, e_end_f1, e_f2, e_f3) in meters

    function Element(nodeStart::Node, nodeEnd::Node, section::Section, id = :element; 
                     release = :fixedfixed, eccentricity = nothing)

        @assert in(release, keys(_ReleaseDict)) "Release not recognized; choose from: :fixedfixed, :freefixed, :fixedfree, :freefree, :joist"

        # Process eccentricity - convert to meters or default to zeros
        if eccentricity === nothing
            ecc = (0.0, 0.0, 0.0, 0.0)
        else
            ecc = (
                to_meters(eccentricity[1]),
                to_meters(eccentricity[2]),
                to_meters(eccentricity[3]),
                to_meters(eccentricity[4])
            )
        end

        element = new{_ReleaseDict[release]}(
            section,
            nodeStart,
            nodeEnd,
            0,
            Vector{Int64}(undef, 12),
            0.0u"m",  # Initialize length as Quantity
            zeros(12, 12),
            zeros(12),
            zeros(12, 12),
            pi/2,
            repeat([zeros(3)], 3),
            zeros(12),
            id,
            ecc
        )

        return element
    end
end

"""
Check if element has eccentric connections.
"""
has_eccentricity(elem::Element) = any(e != 0.0 for e in elem.eccentricity)

"""
    BridgeElement(elementStart, posStart, elementEnd, posEnd, section, id=:element; release=:fixedfixed)

Create a bridge element between two frame elements. 

Connects from `elementStart` at position `elementStart.length * posStart` away from 
`elementStart.nodeStart.position` to `elementEnd` at `elementEnd.length * posEnd` away 
from `elementEnd.nodeStart.position`. 

Note: `posStart, posEnd ∈ ]0, 1[`
"""
mutable struct BridgeElement{R<:Release} <: FrameElement{R}
    elementStart::Element
    posStart::Float64
    elementEnd::Element
    posEnd::Float64
    section::Section
    release::Symbol
    Ψ::Float64
    elementID::Int64
    id::Union{Symbol, Nothing}

    function BridgeElement(elementStart::Element, 
            posStart::Float64, 
            elementEnd::Element, 
            posEnd::Float64, 
            section::Section,
            id = :element;
            release = :fixedfixed)

        @assert 0 < posStart < 1 && 0 < posEnd < 1 "posStart/End must be ∈ ]0,1["
        @assert in(release, keys(_ReleaseDict))

        be = new{_ReleaseDict[release]}(
            elementStart, 
            posStart, 
            elementEnd, 
            posEnd, 
            section, 
            release, 
            pi/2, 
            0, 
            id
        )

        return be
    end
end


"""
    TrussElement(nodeStart, nodeEnd, section, id=:element)

Create a truss element (axial force only, no bending).

# Example
```julia
TrussElement(node1, node2, TrussSection(0.01u"m^2", 200u"GPa"))
```
"""
mutable struct TrussElement <: AbstractElement
    section::Union{TrussSection, Section}   # Cross section
    nodeStart::TrussNode                    # Start node
    nodeEnd::TrussNode                      # End node
    elementID::Int64                        # Element index
    globalID::Vector{Int64}                 # Element global DOFs
    length::QuantityDistance                # Length of element
    K::Matrix{Float64}                      # Stiffness matrix in GCS
    R::Matrix{Float64}                      # Transformation matrix
    forces::Vector{Float64}                 # Elemental forces in LCS
    Ψ::Float64
    LCS::Vector{Vector{Float64}}
    id::Union{Symbol, Nothing}              # Optional identifier

    function TrussElement(nodeStart::TrussNode, nodeEnd::TrussNode, section::AbstractSection, id = :element)
        
        element = new(
            section,
            nodeStart,
            nodeEnd,
            0,
            Vector{Int64}(undef, 6),
            0.0u"m",  # Initialize length as Quantity
            zeros(2, 2),
            zeros(2, 6),
            zeros(6),
            pi/2,
            repeat([zeros(3)], 3),
            id
        )

        return element
    end
end


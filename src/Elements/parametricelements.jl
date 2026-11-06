#=
Parametric Frame Element Definitions for Asap.jl
================================================

`ParametricFrameElement` is a **topological declaration**, not a peer element
type. It records the intent "connect a frame member parametrically between two
host `Element`s at fractional positions along their spans." During `process!`
(via `normalize_parametric_topology!` in `Model/parametricprocessing.jl`) it is resolved
away: intermediate nodes are inserted on the host elements, the hosts are
shattered into sub-`Element`s, and each `ParametricFrameElement` is replaced by
an ordinary `Element` connecting the newly created nodes. After processing, the
model contains only regular `Element`s — no `ParametricFrameElement` survives
into the stiffness assembly.

This type was formerly named `BridgeElement`. The "bridge" name was misleading:
it is not bridge-specific and does not model a distinct structural member. It is
a generic parametric landing of any frame element onto two hosts, which is why
it is now named for what it does (parametric frame connection) rather than a
single application.
=#

"""
    ParametricFrameElement(elementStart, posStart, elementEnd, posEnd, section, id=:element; release=:fixedfixed)

Declare a frame element connected parametrically between two host frame elements.

Connects from `elementStart` at position `elementStart.length * posStart` away from
`elementStart.nodeStart.position` to `elementEnd` at `elementEnd.length * posEnd` away
from `elementEnd.nodeStart.position`.

This is a topological declaration that is resolved into a regular `Element` during
`process!` (see `normalize_parametric_topology!`); it is not a distinct structural
element type.

Note: `posStart, posEnd ∈ ]0, 1[`
"""
mutable struct ParametricFrameElement{R<:Release} <: FrameElement{R}
    elementStart::Element
    posStart::Float64
    elementEnd::Element
    posEnd::Float64
    section::Section
    release::Symbol
    Ψ::Float64
    elementID::Int64
    id::Union{Symbol, Nothing}

    function ParametricFrameElement(elementStart::Element, 
            posStart::Float64, 
            elementEnd::Element, 
            posEnd::Float64, 
            section::Section,
            id = :element;
            release = :fixedfixed)

        @assert 0 < posStart < 1 && 0 < posEnd < 1 "posStart/End must be ∈ ]0,1["
        @assert in(release, keys(_ReleaseDict))

        pfe = new{_ReleaseDict[release]}(
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

        return pfe
    end
end

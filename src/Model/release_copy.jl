# =============================================================================
# Copy of a model with different frame end releases
# =============================================================================

"""Rebuild an element load on `element`, which replaces the load's own element."""
_rebind_load(ld::LineLoad, element::FrameElement) = LineLoad(element, ld.value, ld.id)
_rebind_load(ld::GravityLoad, element::FrameElement) = GravityLoad(element, ld.factor, ld.id)
_rebind_load(ld::PointLoad, element::FrameElement) =
    PointLoad(element, ld.position, ld.value, ld.id)
_rebind_load(ld::TributaryLoad, element::FrameElement) =
    TributaryLoad(element, ld.positions, ld.widths, ld.pressure, ld.direction, ld.id)
_rebind_load(ld::ElementLoad, ::FrameElement) = throw(ArgumentError(
    "with_releases cannot move a $(typeof(ld)) onto a re-released element"))

_set_frame_elements!(model::FrameModel, frames) = (model.elements = frames)
_set_frame_elements!(model::Model, frames) = (model.frame_elements = frames)

"""
    with_releases(model, releases; links = nothing, load_cases = nothing)
        -> (; model, links, load_cases)

An independent copy of `model` in which frame element `i` has the end release
`releases[i]` (`:fixedfixed`, `:fixedfree`, `:freefixed`, `:freefree`, or
`:joist`). Every other element, node, section, support, spring, and load is
copied unchanged. The original model is not modified.

The release is the element's type parameter, so each listed element is rebuilt
with the same nodes, section, `id`, roll `Ψ`, and eccentricity. Element loads on
a rebuilt element (`LineLoad`, `GravityLoad`, `PointLoad`, `TributaryLoad`) are
rebuilt on it with the same values, because a load's fixed-end forces depend on
the release of the element it loads.

`links` (rigid links) and `load_cases` (vectors of loads) refer to the model's
nodes and elements. They are copied with the model, so the returned ones refer
to the copy, and loads in `load_cases` are rebuilt like the model's own loads.
The same original load object maps to the same copied object everywhere.

An `AreaLoad` keeps its shells and pressure, distributes to the rebuilt
elements, and recomputes its tributary split when the copy is processed.

The copy is unprocessed. A model built from parametric frame elements is not
supported.
"""
function with_releases(model::Union{FrameModel, Model},
                       releases::AbstractDict{<:Integer, Symbol};
                       links = nothing, load_cases = nothing)
    model._parametric_source === nothing || throw(ArgumentError(
        "with_releases does not support a model with parametric frame elements"))
    copied, copied_links, copied_cases = deepcopy((model, links, load_cases))
    frames = FrameElement[e for e in _frame_elements(copied)]
    swapped = IdDict{FrameElement, FrameElement}()
    for (i, release) in releases
        1 <= i <= length(frames) || throw(ArgumentError(
            "frame element $i is outside 1:$(length(frames))"))
        haskey(_ReleaseDict, release) || throw(ArgumentError(
            "release :$release is not one of $(sort!(collect(keys(_ReleaseDict))))"))
        old = frames[i]
        old isa Element || throw(ArgumentError(
            "frame element $i is a $(typeof(old)), not an Element"))
        new = Element(old.nodeStart, old.nodeEnd, old.section, old.id; release = release)
        new.Ψ = old.Ψ
        new.eccentricity = old.eccentricity
        frames[i] = new
        swapped[old] = new
    end
    _set_frame_elements!(copied, frames)
    rebuilt = IdDict{AbstractLoad, AbstractLoad}()
    swap(e) = get(swapped, e, e)
    function rebind(ld)
        if ld isa AreaLoad
            # The copy is private, so its beam lists are replaced in place. The
            # tributary split is recomputed from geometry on the next process.
            ld.distribute_to isa Symbol ||
                (ld.distribute_to = FrameElement[swap(e) for e in ld.distribute_to])
            ld.interior_beams = FrameElement[swap(e) for e in ld.interior_beams]
            ld._tributary_loads = nothing
            return ld
        end
        ld isa ElementLoad && haskey(swapped, ld.element) || return ld
        return get!(() -> _rebind_load(ld, swapped[ld.element]), rebuilt, ld)
    end
    copied.loads = AbstractLoad[rebind(ld) for ld in copied.loads]
    if copied_cases !== nothing
        copied_cases = [AbstractLoad[rebind(ld) for ld in case] for case in copied_cases]
    end
    _invalidate_parametric_caches!(copied)
    copied.processed = false
    return (; model = copied, links = copied_links, load_cases = copied_cases)
end

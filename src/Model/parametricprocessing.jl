# =============================================================================
# Parametric frame topology normalization
# =============================================================================

"""
Source topology retained by a model that contains parametric frame declarations.

The source objects are never replaced or processed as analysis elements. Original
nodes are intentionally shared with the realized topology so geometry changes and
mixed-model shell connectivity remain explicit and repeatable.
"""
struct ParametricSource
    nodes::Vector{Node}
    elements::Vector{FrameElement}
    loads::Vector{AbstractLoad}
end

"""Fully validated, off-model realization ready for one atomic commit."""
struct ParametricTopologyResult
    nodes::Vector{Node}
    elements::Vector{FrameElement}
    loads::Vector{AbstractLoad}
    element_map::IdDict{FrameElement,Vector{Element}}
    fracture_fractions::IdDict{Element,Vector{Float64}}
    fracture_nodes::IdDict{Element,Vector{Node}}
end

_frame_elements(model::FrameModel) = model.elements
_frame_elements(model::Model) = model.frame_elements

function _capture_parametric_source(model::Union{FrameModel,Model})
    return ParametricSource(
        copy(model.nodes),
        FrameElement[model_element for model_element in _frame_elements(model)],
        AbstractLoad[load for load in model.loads],
    )
end

_release_symbol(::Element{FixedFixed}) = :fixedfixed
_release_symbol(::Element{FixedFree}) = :fixedfree
_release_symbol(::Element{FreeFixed}) = :freefixed
_release_symbol(::Element{FreeFree}) = :freefree
_release_symbol(::Element{Joist}) = :joist

_start_released(::Element{<:Union{FreeFixed,FreeFree}}) = true
_start_released(::Element) = false
_end_released(::Element{<:Union{FixedFree,FreeFree}}) = true
_end_released(::Element) = false

function _fragment_release(
    host::Element,
    fragment_index::Int,
    fragment_count::Int,
)
    host isa Element{Joist} && throw(ArgumentError(
        "Parametric processing does not support Joist hosts: the current release " *
        "types cannot preserve one-ended bending release with torsional continuity."))
    start_free = fragment_index == 1 && _start_released(host)
    end_free = fragment_index == fragment_count && _end_released(host)
    start_free && end_free && return :freefree
    start_free && return :freefixed
    end_free && return :fixedfree
    return :fixedfixed
end

function _clone_element(
    source::Element,
    node_start::Node,
    node_end::Node;
    release::Symbol = _release_symbol(source),
    first_fragment::Bool = true,
    last_fragment::Bool = true,
)
    eccentricity = (
        first_fragment ? source.eccentricity[1] : 0.0,
        last_fragment ? source.eccentricity[2] : 0.0,
        source.eccentricity[3],
        source.eccentricity[4],
    )
    element = Element(
        node_start, node_end, source.section, source.id;
        release=release,
    )
    element.eccentricity = eccentricity
    element.Ψ = source.Ψ
    return element
end

function _parametric_position(host::Element, fraction::Float64)
    start = Float64[to_meters(value) for value in host.nodeStart.position]
    finish = Float64[to_meters(value) for value in host.nodeEnd.position]
    return (start .+ fraction .* (finish .- start)) .* u"m"
end

function _coalesced_fractions(values::Vector{Float64}, tol::Float64)
    isempty(values) && return Float64[]
    sorted = sort(values)
    unique_values = Float64[sorted[1]]
    for value in @view sorted[2:end]
        abs(value - unique_values[end]) <= tol || push!(unique_values, value)
    end
    return unique_values
end

function _validate_parametric_source(source::ParametricSource)
    declarations = ParametricFrameElement[
        element for element in source.elements
        if element isa ParametricFrameElement
    ]
    isempty(declarations) && throw(ArgumentError(
        "Parametric source contains no ParametricFrameElement declarations."))

    regular = Element[
        element for element in source.elements
        if element isa Element
    ]
    length(regular) + length(declarations) == length(source.elements) ||
        throw(ArgumentError(
            "Parametric processing currently supports only Element and " *
            "ParametricFrameElement frame entries."))

    hosts = IdDict{Element,Bool}(host => true for host in regular)
    for declaration in declarations
        haskey(hosts, declaration.elementStart) || throw(ArgumentError(
            "ParametricFrameElement start host is not present in the source model."))
        haskey(hosts, declaration.elementEnd) || throw(ArgumentError(
            "ParametricFrameElement end host is not present in the source model."))
    end
    return regular, declarations
end

function _clone_node_load(load::NodeForce)
    return NodeForce(load.node, load.value, load.id)
end

function _clone_node_load(load::NodeMoment)
    return NodeMoment(load.node, load.value, load.id)
end

function _remap_line_load(
    load::LineLoad,
    element_map::IdDict{FrameElement,Vector{Element}},
)
    targets = get(element_map, load.element, nothing)
    targets === nothing && throw(ArgumentError(
        "LineLoad references an element outside the parametric source topology."))
    return AbstractLoad[LineLoad(target, load.value, load.id) for target in targets]
end

function _remap_point_load(
    load::PointLoad,
    element_map::IdDict{FrameElement,Vector{Element}},
    fracture_fractions::IdDict{Element,Vector{Float64}},
    fracture_nodes::IdDict{Element,Vector{Node}},
    tol::Float64,
)
    targets = get(element_map, load.element, nothing)
    targets === nothing && throw(ArgumentError(
        "PointLoad references an element outside the parametric source topology."))

    if load.element isa ParametricFrameElement
        return AbstractLoad[PointLoad(only(targets), load.position, load.value, load.id)]
    end

    host = load.element::Element
    fractions = get(fracture_fractions, host, Float64[])
    nodes = get(fracture_nodes, host, Node[])
    for (fraction, node) in zip(fractions, nodes)
        if abs(load.position - fraction) <= tol
            return AbstractLoad[NodeForce(node, load.value, load.id)]
        end
    end

    boundaries = [0.0; fractions; 1.0]
    fragment_index = searchsortedlast(boundaries, load.position)
    fragment_index = min(fragment_index, length(targets))
    lower = boundaries[fragment_index]
    upper = boundaries[fragment_index + 1]
    local_fraction = (load.position - lower) / (upper - lower)
    return AbstractLoad[
        PointLoad(targets[fragment_index], local_fraction, load.value, load.id)
    ]
end

function _remap_loads(
    source::ParametricSource,
    element_map::IdDict{FrameElement,Vector{Element}},
    fracture_fractions::IdDict{Element,Vector{Float64}},
    fracture_nodes::IdDict{Element,Vector{Node}},
    tol::Float64,
)
    remapped = AbstractLoad[]
    for load in source.loads
        if load isa NodeForce || load isa NodeMoment
            push!(remapped, _clone_node_load(load))
        elseif load isa LineLoad
            append!(remapped, _remap_line_load(load, element_map))
        elseif load isa PointLoad
            append!(remapped, _remap_point_load(
                load, element_map, fracture_fractions, fracture_nodes, tol))
        else
            throw(ArgumentError(
                "Parametric processing does not support $(nameof(typeof(load))); " *
                "supported loads are NodeForce, NodeMoment, LineLoad, and PointLoad."))
        end
    end
    return remapped
end

"""
    realize_parametric_topology(source, tol) -> ParametricTopologyResult

Build and validate a complete realized topology without mutating a model.
"""
function realize_parametric_topology(source::ParametricSource, tol::Float64)
    regular, declarations = _validate_parametric_source(source)

    requested = IdDict{Element,Vector{Float64}}()
    for declaration in declarations
        push!(get!(requested, declaration.elementStart, Float64[]), declaration.posStart)
        push!(get!(requested, declaration.elementEnd, Float64[]), declaration.posEnd)
    end

    fracture_fractions = IdDict{Element,Vector{Float64}}()
    fracture_nodes = IdDict{Element,Vector{Node}}()
    element_map = IdDict{FrameElement,Vector{Element}}()
    realized_nodes = copy(source.nodes)
    realized_elements = FrameElement[]

    for host in regular
        fractions = _coalesced_fractions(get(requested, host, Float64[]), tol)
        fracture_fractions[host] = fractions
        nodes = Node[]
        for fraction in fractions
            node = Node(_parametric_position(host, fraction), :free)
            node.id = :parametric
            push!(nodes, node)
            push!(realized_nodes, node)
        end
        fracture_nodes[host] = nodes

        chain = Node[host.nodeStart; nodes; host.nodeEnd]
        fragments = Element[]
        fragment_count = length(chain) - 1
        for index in 1:fragment_count
            fragment = _clone_element(
                host, chain[index], chain[index + 1];
                release=_fragment_release(host, index, fragment_count),
                first_fragment=index == 1,
                last_fragment=index == fragment_count,
            )
            push!(fragments, fragment)
            push!(realized_elements, fragment)
        end
        element_map[host] = fragments
    end

    for declaration in declarations
        start_fractions = fracture_fractions[declaration.elementStart]
        end_fractions = fracture_fractions[declaration.elementEnd]
        start_index = argmin(abs.(start_fractions .- declaration.posStart))
        end_index = argmin(abs.(end_fractions .- declaration.posEnd))
        bridge = Element(
            fracture_nodes[declaration.elementStart][start_index],
            fracture_nodes[declaration.elementEnd][end_index],
            declaration.section,
            declaration.id;
            release=declaration.release,
        )
        bridge.Ψ = declaration.Ψ
        push!(realized_elements, bridge)
        element_map[declaration] = Element[bridge]
    end

    loads = _remap_loads(
        source, element_map, fracture_fractions, fracture_nodes, tol)
    return ParametricTopologyResult(
        realized_nodes,
        realized_elements,
        loads,
        element_map,
        fracture_fractions,
        fracture_nodes,
    )
end

function _validate_spring_nodes(model, result::ParametricTopologyResult)
    realized = IdDict{Node,Bool}(node => true for node in result.nodes)
    for spring in model.springs
        haskey(realized, spring.node) || throw(ArgumentError(
            "A model-owned spring references a generated or removed node. " *
            "Attach persistent springs to source nodes."))
    end
end

function _commit_parametric_topology!(
    model::FrameModel,
    result::ParametricTopologyResult,
)
    model.nodes = result.nodes
    model.elements = result.elements
    model.loads = result.loads
    model.nNodes = length(result.nodes)
    model.nElements = length(result.elements)
    _invalidate_parametric_caches!(model)
end

function _commit_parametric_topology!(
    model::Model,
    result::ParametricTopologyResult,
)
    model.nodes = result.nodes
    model.frame_elements = result.elements
    model.loads = result.loads
    model.nNodes = length(result.nodes)
    model.nFrameElements = length(result.elements)
    model.nShellElements = length(model.shell_elements)
    _invalidate_parametric_caches!(model)
end

function _invalidate_parametric_caches!(model::Union{FrameModel,Model})
    model._factorization = nothing
    model._factorization_stale = false
    model._factor_key = UInt(0)
    model._scatter = nothing
    model._elemental_loads = nothing
end

"""
    normalize_parametric_topology!(model) -> Bool

Atomically realize retained parametric source topology. Returns `true` when the
model owns a parametric source and was normalized.

Supported remapping is deliberately limited to `NodeForce`, `NodeMoment`,
`LineLoad`, and `PointLoad`. `TributaryLoad` remapping, joist-host fracture,
and coordinate-based frame/shell node welding are deferred until each has an
explicit conservation or release-preservation contract.
"""
function normalize_parametric_topology!(model::Union{FrameModel,Model})
    source = if model._parametric_source === nothing
        any(element -> element isa ParametricFrameElement, _frame_elements(model)) ||
            return false
        _capture_parametric_source(model)
    else
        model._parametric_source::ParametricSource
    end

    result = realize_parametric_topology(source, model.tol)
    _validate_spring_nodes(model, result)
    _commit_parametric_topology!(model, result)
    model._parametric_source = source
    return true
end


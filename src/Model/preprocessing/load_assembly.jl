# =============================================================================
# Load Population
# =============================================================================

"""
    populate_load!(model::AbstractModel, load::NodeForce)

Populate the global load vector `model.P` with a nodal force.
"""
function populate_load!(model::AbstractModel, load::NodeForce)
    gid = load.node.globalID
    @inbounds for i in 1:3
        model.P[gid[i]] += to_newtons(load.value[i])
    end
end

"""
    populate_load!(P::Vector{Float64}, load::NodeForce)

Populate the global load vector P with a nodal force.
"""
function populate_load!(P::Vector{Float64}, load::NodeForce)
    gid = load.node.globalID
    @inbounds for i in 1:3
        P[gid[i]] += to_newtons(load.value[i])
    end
end

"""
    populate_load!(model, load::NodeMoment)

Populate the global load vector `model.P` with a nodal moment.
"""
function populate_load!(model::AbstractModel, load::NodeMoment)
    gid = load.node.globalID
    @inbounds for i in 1:3
        model.P[gid[i+3]] += to_newton_meters(load.value[i])
    end
end

"""
    populate_load!(P::Vector{Float64}, load::NodeMoment)

Populate the global load vector P with a nodal moment.
"""
function populate_load!(P::Vector{Float64}, load::NodeMoment)
    gid = load.node.globalID
    @inbounds for i in 1:3
        P[gid[i+3]] += to_newton_meters(load.value[i])
    end
end

"""
    populate_load!(model::FrameModel, load::ElementLoad)

Generate the fixed-end force vector `Q` for a given load, and populate the global fixed-end force vector `Pf`.
"""
function populate_load!(model::FrameModel, load::ElementLoad)
    idx = load.element.globalID
    q_local = q(load)
    buf = _thread_buffer(_POPULATE_LOAD_BUF)
    mul!(buf, transpose(load.element.R), q_local)
    @inbounds for i in eachindex(idx)
        load.element.Q[i] += buf[i]
        model.Pf[idx[i]] += buf[i]
    end
end

"""
    populate_load!(model::Model, load::ElementLoad)

Generate the fixed-end force vector for a given load in a unified model.
"""
function populate_load!(model::Model, load::ElementLoad)
    idx = load.element.globalID
    q_local = q(load)
    buf = _thread_buffer(_POPULATE_LOAD_BUF)
    mul!(buf, transpose(load.element.R), q_local)
    @inbounds for i in eachindex(idx)
        load.element.Q[i] += buf[i]
        model.Pf[idx[i]] += buf[i]
    end
end

"""
    populate_load!(Pf::Vector{Float64}, load::ElementLoad)

Populate an external fixed-end force vector Pf with respect to an elemental load.
"""
function populate_load!(Pf::Vector{Float64}, load::ElementLoad)
    idx = load.element.globalID
    q_local = q(load)
    buf = _thread_buffer(_POPULATE_LOAD_BUF)
    mul!(buf, transpose(load.element.R), q_local)
    @inbounds for i in eachindex(idx)
        Pf[idx[i]] += buf[i]
    end
end

"""
    populate_load!(model, load::AreaLoad)

Populate the global load vector with forces from an area load.

Behavior depends on `load.distribute_to`:
- `:nodes` → FEM approach: equivalent nodal forces on shell nodes
- `Vector{FrameElement}` → Tributary approach: compute tributary loads for beams
"""
function populate_load!(model::AbstractModel, load::AreaLoad)
    if load.distribute_to == :nodes
        # FEM approach: inline nodal forces — zero allocations
        p = ustrip(u"Pa", load.pressure)
        d1, d2, d3 = load.direction
        _accumulate_shell_pressure!(model.P, load.shells, p, d1, d2, d3)
    else
        # Tributary approach: geometry cached, only pressure updated
        if isnothing(load._tributary_loads)
            axis_vec = isnothing(load.axis) ? nothing : [load.axis[1], load.axis[2]]
            load._tributary_loads = _shell_to_tributary_loads(
                load.shells, 
                load.distribute_to, 
                load.pressure;
                axis=axis_vec,
                direction=load.direction,
                interior_beams=load.interior_beams
            )
        else
            # Geometry unchanged — just update pressure on cached loads
            for trib_load in load._tributary_loads
                trib_load.pressure = load.pressure
            end
        end
        
        # Apply each tributary load
        for trib_load in load._tributary_loads
            populate_load!(model, trib_load)
        end
    end
end

"""
Add an area load's nodal forces. `shells` is stored as
`Vector{<:ShellElement}`; the barrier makes `area` and `nodes` concrete.
"""
function _accumulate_shell_pressure!(
    P::Vector{Float64}, shells::Vector{S}, p::Float64, d1::Float64, d2::Float64, d3::Float64,
) where {S}
    @inbounds for shell in shells
        fpn = p * shell.area / length(shell.nodes)
        for node in shell.nodes
            gid = node.globalID
            P[gid[1]] += fpn * d1
            P[gid[2]] += fpn * d2
            P[gid[3]] += fpn * d3
        end
    end
    return P
end

"""
    populate_load!(P::Vector{Float64}, load::AreaLoad)

Populate an external load vector with area load contributions.
"""
function populate_load!(P::Vector{Float64}, load::AreaLoad)
    if load.distribute_to == :nodes
        for (node, fvec) in nodal_forces(load)
            gid = node.globalID
            @inbounds for k in 1:3
                P[gid[k]] += fvec[k]
            end
        end
    else
        error("Tributary distribution requires model context. Use populate_load!(model, load) instead.")
    end
end

# =============================================================================
# Load Vector Assembly
# =============================================================================

"""
    populate_loads!(model::FrameModel)

Generate the nodal force vectors `model.P` (external) and `model.Pf` (fixed-end)
"""
function populate_loads!(model::FrameModel)
    fill!(model.P, 0.0)
    fill!(model.Pf, 0.0)

    for load in model.loads
        populate_load!(model, load)
    end
end

"""
    populate_loads!(model::ShellModel)

Generate the nodal force vectors `model.P` for a shell model.
"""
function populate_loads!(model::ShellModel)
    fill!(model.P, 0.0)

    for load in model.loads
        populate_load!(model, load)
    end
end

"""
    populate_loads!(model::Model)

Generate the nodal force vectors for a unified model.
"""
function populate_loads!(model::Model)
    fill!(model.P, 0.0)
    fill!(model.Pf, 0.0)

    for load in model.loads
        populate_load!(model, load)
    end
end

"""
    populate_loads!(model::TrussModel)

Generate the nodal force vectors `model.P`
"""
function populate_loads!(model::TrussModel)
    fill!(model.P, 0.0)

    for load in model.loads
        populate_load!(model, load)
    end
end


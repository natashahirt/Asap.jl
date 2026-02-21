"""
    planarize!(model::AbstractModel, plane = :XY)

Fix all nodal DOFs to remain on plane = `plane`
"""
function planarize!(model::FrameModel, plane = :XY)
    planarize!(model.nodes, plane)
    if plane == :XY
        for element in model.elements
            element.Ψ = 0.
        end
    end
end

function planarize!(model::Model, plane = :XY)
    planarize!(model.nodes, plane)
    if plane == :XY
        for element in model.frame_elements
            element.Ψ = 0.
        end
    end
end

function planarize!(model::ShellModel, plane = :XY)
    planarize!(model.nodes, plane)
end

function planarize!(model::TrussModel, plane = :XY)
    planarize!(model.nodes, plane)
end

"""
    connectivity(model::AbstractModel)

Get the [nₑ × nₙ] sparse matrix where C[i, j] = -1 if element i starts at node j, 
and C[i,j] = 1 if element i ends at node j, and 0 otherwise.

Only works for frame/truss models with 2-node elements.
"""
function connectivity(model::FrameModel)
    I = vcat([[i, i] for i = 1:model.nElements]...)
    J = vcat([nodeids(e) for e in model.elements]...)
    V = repeat([-1, 1], model.nElements)
    return sparse(I, J, V)
end

function connectivity(model::Model)
    if isempty(model.frame_elements)
        error("connectivity() only works for frame elements")
    end
    n_elem = model.nFrameElements
    I = vcat([[i, i] for i = 1:n_elem]...)
    J = vcat([nodeids(e) for e in model.frame_elements]...)
    V = repeat([-1, 1], n_elem)
    return sparse(I, J, V)
end

function connectivity(model::TrussModel)
    I = vcat([[i, i] for i = 1:model.nElements]...)
    J = vcat([nodeids(e) for e in model.elements]...)
    V = repeat([-1, 1], model.nElements)
    return sparse(I, J, V)
end

"""
    node_positions(model::AbstractModel)

Generate the [nₙ × 3] node position matrix
"""
function node_positions(model::AbstractModel)
    return vcat([[to_meters(p) for p in node.position]' for node in model.nodes]...)
end

export nodePositions
nodePositions(model::AbstractModel) = node_positions(model)

"""
    update_DOF!(model::AbstractModel)

Update the free/fixed degrees of freedom for a model
"""
function update_DOF!(model::AbstractModel)
    n_dof_per_node = model isa TrussModel ? 3 : 6
    _populate_DOFs_flat!(model, length(model.nodes), n_dof_per_node)
end

export updateDOF!
updateDOF!(model::AbstractModel) = update_DOF!(model)

"""
    volume(model::AbstractModel)

Get the material volume of a structural model.
"""
function volume(model::FrameModel)
    lengths = [to_meters(el.length) for el in model.elements]
    areas = [to_meters_squared(sec.A) for sec in getproperty.(model.elements, :section)]
    dot(lengths, areas)
end

function volume(model::Model)
    vol = 0.0
    
    # Frame elements
    if !isempty(model.frame_elements)
        lengths = [to_meters(el.length) for el in model.frame_elements]
        areas = [to_meters_squared(sec.A) for sec in getproperty.(model.frame_elements, :section)]
        vol += dot(lengths, areas)
    end
    
    # Shell elements
    for elem in model.shell_elements
        vol += elem.area * elem.thickness
    end
    
    return vol
end

function volume(model::ShellModel)
    vol = 0.0
    for elem in model.elements
        vol += elem.area * elem.thickness
    end
    return vol
end

function volume(model::TrussModel)
    lengths = [to_meters(el.length) for el in model.elements]
    areas = [to_meters_squared(sec.A) for sec in getproperty.(model.elements, :section)]
    dot(lengths, areas)
end

function Base.copy(model::TrussModel)
    
    #new model
    nodes = Vector{TrussNode}()
    elements = Vector{TrussElement}()
    loads = Vector{NodeForce}()

    #new nodes
    for node in model.nodes
        newnode = TrussNode(deepcopy(node.position), node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for element in model.elements
        newelement = TrussElement(nodes, deepcopy(element.nodeIDs), element.section)
        newelement.id = element.id
        push!(elements, newelement)
    end

    #new loads
    for load in model.loads
        newload = NodeForce(nodes[load.node.nodeID], deepcopy(load.value))
        newload.id = load.id
        push!(loads, newload)
    end

    model = TrussModel(nodes, elements, loads)
    solve!(model)

    return model
end

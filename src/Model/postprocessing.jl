"""
    reactions!(model::AbstractModel)

Populate external reaction forces in `model.reactions`
"""
function reactions!(model::Model)
    model.reactions = zeros(model.nDOFs)
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u + model.Pf[model.fixedDOFs]
end

function reactions!(model::TrussModel)
    model.reactions = zeros(model.nDOFs)
    model.reactions[model.fixedDOFs] = model.S[model.fixedDOFs, :] * model.u
end

"""
    post_process_nodes!(model::AbstractModel)

Populate nodal reaction and displacement fields with Unitful quantities.
Reactions are in SI base units: forces in N, moments in N*m.
Displacements are in SI base units: translations in m, rotations in radians.
"""
function post_process_nodes!(model::Model)
    for node in model.nodes
        global_ids = node.globalID
        # Convert reactions: DOFs 1-3 are forces (N), DOFs 4-6 are moments (N*m)
        # Convert displacements: DOFs 1-3 are translations (m), DOFs 4-6 are rotations (rad)
        reactions_quantities = Vector{Quantity}(undef, length(global_ids))
        displacements_quantities = Vector{Quantity}(undef, length(global_ids))
        
        for (i, gid) in enumerate(global_ids)
            # Reactions: forces for DOFs 1-3, moments for DOFs 4-6
            if i <= 3
                reactions_quantities[i] = model.reactions[gid] * u"N"
            else
                reactions_quantities[i] = model.reactions[gid] * u"N*m"
            end
            
            # Displacements: translations for DOFs 1-3, rotations for DOFs 4-6
            if i <= 3
                displacements_quantities[i] = model.u[gid] * u"m"
            else
                displacements_quantities[i] = model.u[gid] * u"rad"
            end
        end
        
        node.reaction = reactions_quantities
        node.displacement = displacements_quantities
    end
end

function post_process_nodes!(model::TrussModel)
    for node in model.nodes
        global_ids = node.globalID
        # Truss: all 3 DOFs are forces (N) and translations (m)
        reactions_quantities = [model.reactions[gid] * u"N" for gid in global_ids]
        displacements_quantities = [model.u[gid] * u"m" for gid in global_ids]
        
        node.reaction = reactions_quantities
        node.displacement = displacements_quantities
    end
end

"""
    post_process_elements!(model::AbstractModel)

Populate elemental LCS force vectors in `element.forces`
"""
function post_process_elements!(model::Model)
    for element in model.elements
        element.forces = element.R * (element.K * model.u[element.globalID] + element.Q)
    end
end

function post_process_elements!(model::TrussModel)
    for element in model.elements
        element.forces = element.R * (element.K * model.u[element.globalID])
    end
end

"""
    post_process!(model::AbstractModel)

Post process a model after solving for displacements using `solve!(model)`
"""
function post_process!(model::AbstractModel)
    reactions!(model)
    post_process_nodes!(model)
    post_process_elements!(model)
end
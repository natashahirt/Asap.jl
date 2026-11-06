# =============================================================================
# Load Vector Creation (for external load vectors)
# =============================================================================

"""
    create_F(model::FrameModel, loads::Vector{AbstractLoad})

Create load vector F = P - Pf for a given load set.
"""
function create_F(model::FrameModel, loads::Vector{T}) where {T<:AbstractLoad}
    P = zeros(model.nDOFs)
    Pf = zeros(model.nDOFs)

    for load in loads
        if load isa ElementLoad
            populate_load!(Pf, load)
        else
            populate_load!(P, load)
        end
    end

    return P - Pf
end

"""
    create_F(model::Model, loads::Vector{AbstractLoad})

Create load vector F = P - Pf for a given load set.
"""
function create_F(model::Model, loads::Vector{T}) where {T<:AbstractLoad}
    P = zeros(model.nDOFs)
    Pf = zeros(model.nDOFs)

    for load in loads
        if load isa ElementLoad
            populate_load!(Pf, load)
        else
            populate_load!(P, load)
        end
    end

    return P - Pf
end

"""
    create_F(model::TrussModel, loads::Vector{NodeForce})

Create load vector F = P for truss models.
"""
function create_F(model::TrussModel, loads::Vector{NodeForce})
    P = zeros(model.nDOFs)

    for load in loads
        populate_load!(P, load)
    end

    return P
end


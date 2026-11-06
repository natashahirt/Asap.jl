# =============================================================================
# Solving
# =============================================================================

"""
    solve!(model::FrameModel; reprocess = false)

Solve for the nodal displacements of a structural model. 
`reprocess = true` reevaluates all node/element properties and reassembles the global stiffness matrix.
"""
function solve!(model::FrameModel; reprocess = false, postprocess::Symbol = :all)
    if !model.processed || reprocess
        for element in model.elements
            if element isa Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    fact = _get_factorization(model)
    U = fact \ F

    model.compliance = dot(U, F)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

"""
    solve!(model::ShellModel; reprocess = false)

Solve for the nodal displacements of a shell model.
"""
function solve!(model::ShellModel; reprocess = false, postprocess::Symbol = :all)
    if !model.processed || reprocess
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    # CHOLMOD's Julia 1.12 solve path does not accept a non-strided indexed
    # `SubArray`; materialize the free-load vector at the solver boundary.
    P_free = model.P[idx]
    fact = _get_factorization(model)
    U = fact \ P_free

    model.compliance = dot(U, P_free)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

"""
    solve!(model::Model; reprocess = false)

Solve for the nodal displacements of a unified model.
Works with mixed frame+shell models or single-type models.
"""
function solve!(model::Model; reprocess = false, postprocess::Symbol = :all)
    if !model.processed || reprocess
        for element in model.frame_elements
            if element isa Element
                element.Q = zero(element.Q)
            end
        end
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    F = model.P[idx] - model.Pf[idx]
    fact = _get_factorization(model)
    U = fact \ F

    model.compliance = dot(U, F)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end

"""
    solve!(model::TrussModel; reprocess = false)

Solve for the nodal displacements of a structural truss model.
"""
function solve!(model::TrussModel; reprocess = false, postprocess::Symbol = :all)
    if !model.processed || reprocess
        process!(model)
        model._factorization = nothing
    end

    idx = model.freeDOFs
    # CHOLMOD's Julia 1.12 solve path requires a materialized vector.
    P_free = model.P[idx]
    fact = _get_factorization(model)
    U = fact \ P_free

    model.compliance = dot(U, P_free)
    if length(model.u) == model.nDOFs; fill!(model.u, 0.0) else model.u = zeros(model.nDOFs) end
    model.u[idx] = U

    post_process!(model; targets=postprocess)
end


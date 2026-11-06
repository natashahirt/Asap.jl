# =============================================================================
# Solving with Custom Loads
# =============================================================================

"""
    solve(model::FrameModel, L::Vector{AbstractLoad})

Return the displacement vector under a given load set `L`, reusing the cached
stiffness factorization.
"""
function solve(model::FrameModel, L::Vector{<:AbstractLoad})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    fact = _get_factorization(model)
    
    u = zeros(model.nDOFs)
    u[idx] = fact \ F[idx]
    return u
end

"""
    solve(model::Model, L::Vector{AbstractLoad})

Return the displacement vector under a given load set `L`, reusing the cached
stiffness factorization.
"""
function solve(model::Model, L::Vector{<:AbstractLoad})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    fact = _get_factorization(model)
    
    u = zeros(model.nDOFs)
    u[idx] = fact \ F[idx]
    return u
end

"""
    solve(model::FrameModel, cases::Vector{<:Vector{<:AbstractLoad}})

Solve multiple load cases against the same stiffness matrix. Returns a vector
of displacement vectors, one per case. Uses a single factorization and a
multi-RHS LAPACK solve for efficiency.
"""
function solve(model::FrameModel, cases::Vector{<:Vector{<:AbstractLoad}})
    model.processed || process!(model)
    fact = _get_factorization(model)
    idx = model.freeDOFs
    n = model.nDOFs
    ncases = length(cases)

    # Build force matrix — each column is a load case
    F = Matrix{Float64}(undef, length(idx), ncases)
    for (j, loads) in enumerate(cases)
        Fj = create_F(model, loads)
        F[:, j] = Fj[idx]
    end

    # Multi-RHS solve (LAPACK level-3 BLAS)
    U = fact \ F

    # Unpack into full DOF vectors
    results = Vector{Vector{Float64}}(undef, ncases)
    for j in 1:ncases
        u = zeros(n)
        u[idx] = @view U[:, j]
        results[j] = u
    end
    return results
end

"""
    solve(model::Model, cases::Vector{<:Vector{<:AbstractLoad}})

Solve multiple load cases against the same stiffness matrix. Returns a vector
of displacement vectors, one per case. Uses a single factorization and a
multi-RHS LAPACK solve for efficiency.
"""
function solve(model::Model, cases::Vector{<:Vector{<:AbstractLoad}})
    model.processed || process!(model)
    fact = _get_factorization(model)
    idx = model.freeDOFs
    n = model.nDOFs
    ncases = length(cases)

    # Build force matrix — each column is a load case
    F = Matrix{Float64}(undef, length(idx), ncases)
    for (j, loads) in enumerate(cases)
        Fj = create_F(model, loads)
        F[:, j] = Fj[idx]
    end

    # Multi-RHS solve (LAPACK level-3 BLAS)
    U = fact \ F

    # Unpack into full DOF vectors
    results = Vector{Vector{Float64}}(undef, ncases)
    for j in 1:ncases
        u = zeros(n)
        u[idx] = @view U[:, j]
        results[j] = u
    end
    return results
end

"""
    solve!(model::FrameModel, L::Vector{AbstractLoad})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::FrameModel, L::Vector{<:AbstractLoad})
    model.loads = L
    model._elemental_loads = nothing
    process!(model)
    solve!(model)
    post_process!(model)
end

"""
    solve!(model::Model, L::Vector{AbstractLoad})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::Model, L::Vector{<:AbstractLoad})
    model.loads = L
    model._elemental_loads = nothing
    process!(model)
    solve!(model)
    post_process!(model)
end

"""
    solve(model::TrussModel, L::Vector{NodeForce})

Return the displacement vector under a given load set `L`, reusing the cached
stiffness factorization.
"""
function solve(model::TrussModel, L::Vector{NodeForce})
    model.processed || process!(model)
    
    F = create_F(model, L)
    idx = model.freeDOFs
    fact = _get_factorization(model)
    
    u = zeros(model.nDOFs)
    u[idx] = fact \ F[idx]
    return u
end

"""
    solve!(model::TrussModel, L::Vector{NodeForce})

Replace the assigned model loads with a new load vector and solve.
"""
function solve!(model::TrussModel, L::Vector{NodeForce})
    model.loads = L
    process!(model)
    solve!(model)
    post_process!(model)
end


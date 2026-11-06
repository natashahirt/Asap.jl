# =============================================================================
# TrussModel - unchanged
# =============================================================================

"""
    TrussModel(nodes::Vector{TrussNode}, elements::Vector{TrussElement}, loads::Vector{NodeForce})

Create a complete structural model ready for analysis.
"""
mutable struct TrussModel <: ElementModel
    nodes::Vector{TrussNode}
    elements::Vector{TrussElement}
    loads::Vector{NodeForce}
    nNodes::Int64
    nElements::Int64
    DOFs::Vector{Bool}
    nDOFs::Int64
    freeDOFs::Vector{Int64}
    fixedDOFs::Vector{Int64}
    S::SparseMatrixCSC{Float64,Int64}
    P::Vector{Float64}
    u::Vector{Float64}
    reactions::Vector{Float64}
    compliance::Float64
    tol::Float64
    processed::Bool
    _factorization::Union{Nothing, Factorization}  # Cached stiffness factorization
    _scatter::Union{Nothing, _ScatterMap}  # Cached nzval positions for in-place S updates
    _factorization_stale::Bool  # S values changed in place; refactor reusing the symbolic analysis
    _factor_key::UInt           # Sparsity-pattern key of `_factorization` (0 = not reusable)
    
    function TrussModel(nodes::Vector{TrussNode}, elements::Vector{TrussElement}, loads::Vector{NodeForce})
        nnodes = length(nodes)
        nelements = length(elements)

        structure = new(
            nodes,
            elements,
            loads,
            nnodes,
            nelements,
            Bool[],
            0,
            Int64[],
            Int64[],
            spzeros(Float64, 3nnodes, 3nnodes),
            zeros(3nnodes),
            zeros(3nnodes),
            zeros(3nnodes),
            0.0,
            1e-6,
            false,
            nothing,
            nothing,
            false,
            UInt(0)
        )

        return structure
    end
end


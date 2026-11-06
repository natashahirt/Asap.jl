# =============================================================================
# ShellModel - For shell elements only (slabs, walls)
# =============================================================================

"""
    ShellModel(nodes, elements, loads)

Create a structural model with only shell elements (slabs/walls).

For mixed frame+shell models, use `Model` instead.
"""
mutable struct ShellModel{E<:ShellElement,L<:AbstractLoad} <: ElementModel
    nodes::Vector{Node}
    elements::Vector{E}
    loads::Vector{L}
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
    
    function ShellModel(nodes::Vector{Node}, elements::Vector{E}, loads::Vector{L}) where {E<:ShellElement, L<:AbstractLoad}
        nnodes = length(nodes)
        nelements = length(elements)

        structure = new{E,L}(
            nodes,
            elements,
            loads,
            nnodes,
            nelements,
            Bool[],
            0,
            Int64[],
            Int64[],
            spzeros(Float64, 6nnodes, 6nnodes),
            zeros(6nnodes),
            zeros(6nnodes),
            zeros(6nnodes),
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


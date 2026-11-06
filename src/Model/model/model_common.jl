abstract type AbstractModel end
abstract type ElementModel <: AbstractModel end  # Models with single .elements field

"""
    _ScatterMap

Position in `S.nzval` of every element-stiffness entry, in the COO order used by
`create_S!` (elements in assembly order, then `i` outer / `j` inner). Lets
`_update_S_values!` re-assemble `S` with direct array adds instead of a sparse
binary search per entry. Valid only while `model.S === S` and `nnz(S) == nnz`.
"""
struct _ScatterMap
    S::SparseMatrixCSC{Float64,Int64}
    nnz::Int
    idx::Vector{Int}
end

# =============================================================================
# Shared Model Initialization
# =============================================================================

"""Default model state values for initializing any model type."""
function _default_model_state(n_dof::Int)
    (
        DOFs = Bool[],
        nDOFs = 0,
        freeDOFs = Int64[],
        fixedDOFs = Int64[],
        S = spzeros(Float64, n_dof, n_dof),
        P = zeros(n_dof),
        u = zeros(n_dof),
        reactions = zeros(n_dof),
        compliance = 0.0,
        tol = 1e-6,
        processed = false
    )
end

# =============================================================================
# make_ids! - Assign sequential IDs to model components
# =============================================================================

function make_ids!(nodes::Vector{<:AbstractNode})
    for (i, node) in enumerate(nodes)
        node.nodeID = i
    end
end

function make_ids!(elements::Vector{TrussElement})
    for (i, element) in enumerate(elements)
        element.elementID = i
    end
end

function make_ids!(elements::Vector{<:FrameElement})
    for (i, element) in enumerate(elements)
        element.elementID = i
    end
end

function make_ids!(elements::Vector{<:ShellElement})
    for (i, element) in enumerate(elements)
        element.elementID = i
    end
end

function make_ids!(loads::Vector{<:AbstractLoad})
    for (i, load) in enumerate(loads)
        load.loadID = i
    end
end


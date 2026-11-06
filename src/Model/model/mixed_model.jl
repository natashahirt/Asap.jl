# =============================================================================
# Model - Unified model for mixed frame + shell elements
# =============================================================================

"""
    Model(nodes, frame_elements, shell_elements, loads)
    Model(nodes, frame_elements, loads)  # Frame-only (dispatches to FrameModel internally)
    Model(nodes, shell_elements, loads)  # Shell-only (dispatches to ShellModel internally)

Create a unified structural model that can contain both frame elements 
(beams, columns) and shell elements (slabs, walls).

# Examples
```julia
# Mixed model with frames and shells
model = Model(nodes, beams, slabs, loads)

# Frame-only model
model = Model(nodes, beams, loads)

# Shell-only model  
model = Model(nodes, slabs, loads)
```
"""
mutable struct Model <: AbstractModel
    nodes::Vector{Node}
    frame_elements::Vector{<:FrameElement}
    shell_elements::Vector{<:ShellElement}
    loads::Vector{<:AbstractLoad}
    nNodes::Int64
    nFrameElements::Int64
    nShellElements::Int64
    DOFs::Vector{Bool}
    nDOFs::Int64
    freeDOFs::Vector{Int64}
    fixedDOFs::Vector{Int64}
    S::SparseMatrixCSC{Float64,Int64}
    P::Vector{Float64}
    Pf::Vector{Float64}
    u::Vector{Float64}
    reactions::Vector{Float64}
    compliance::Float64
    tol::Float64
    processed::Bool
    _factorization::Union{Nothing, Factorization}  # Cached stiffness factorization
    _elemental_loads::Union{Nothing, Vector{Vector{AbstractLoad}}}  # Lazy cache
    _scatter::Union{Nothing, _ScatterMap}  # Cached nzval positions for in-place S updates
    _factorization_stale::Bool  # S values changed in place; refactor reusing the symbolic analysis
    _factor_key::UInt           # Sparsity-pattern key of `_factorization` (0 = not reusable)
    _parametric_source::Any      # Immutable source-topology references, realized before assembly
    springs::Vector{Spring}      # Model-owned grounded springs, reassembled after every process!
    
    # Full constructor: frame + shell elements
    function Model(
        nodes::Vector{Node},
        frame_elements::Vector{E1},
        shell_elements::Vector{E2},
        loads::Vector{L}
    ) where {E1<:FrameElement, E2<:ShellElement, L<:AbstractLoad}
        
        # Merge shell nodes into nodes list (shells may have interior mesh nodes)
        all_nodes = _merge_shell_nodes(nodes, shell_elements)
        
        nnodes = length(all_nodes)
        n_frame = length(frame_elements)
        n_shell = length(shell_elements)

        new(
            all_nodes,
            frame_elements,
            shell_elements,
            loads,
            nnodes,
            n_frame,
            n_shell,
            Bool[],
            0,
            Int64[],
            Int64[],
            spzeros(Float64, 6nnodes, 6nnodes),
            zeros(6nnodes),
            zeros(6nnodes),
            zeros(6nnodes),
            zeros(6nnodes),
            0.0,
            1e-6,
            false,
            nothing,
            nothing,
            nothing,
            false,
            UInt(0),
            nothing,
            Spring[]
        )
    end
end

"""
    _merge_shell_nodes(nodes, shell_elements)

Merge shell element nodes into the nodes list, adding any interior mesh nodes
that aren't already present. Uses object identity to avoid duplicates.
"""
function _merge_shell_nodes(nodes::Vector{Node}, shell_elements::Vector{<:ShellElement})
    isempty(shell_elements) && return nodes
    
    # Track existing nodes by identity
    seen = Set{UInt64}(objectid(n) for n in nodes)
    merged = copy(nodes)
    
    # Add unique shell nodes not already in the list
    for elem in shell_elements
        for node in elem.nodes
            id = objectid(node)
            if !(id in seen)
                push!(seen, id)
                push!(merged, node)
            end
        end
    end
    
    return merged
end

# Convenience constructor: Frame elements only
function Model(
    nodes::Vector{Node},
    frame_elements::Vector{E},
    loads::Vector{L}
) where {E<:FrameElement, L<:AbstractLoad}
    # Create Model with empty shell elements
    Model(nodes, frame_elements, ShellTri3[], loads)
end

# Convenience constructor: Shell elements only
function Model(
    nodes::Vector{Node},
    shell_elements::Vector{E},
    loads::Vector{L}
) where {E<:ShellElement, L<:AbstractLoad}
    # Create Model with empty frame elements
    Model(nodes, Element{FixedFixed}[], shell_elements, loads)
end

# =============================================================================
# Helper functions for Model
# =============================================================================

"""Check if model has frame elements."""
has_frame_elements(model::Model) = !isempty(model.frame_elements)

"""Check if model has shell elements."""
has_shell_elements(model::Model) = !isempty(model.shell_elements)

"""Check if model is mixed (has both frame and shell elements)."""
is_mixed(model::Model) = has_frame_elements(model) && has_shell_elements(model)

"""Get all elements (frame + shell) as a vector."""
function all_elements(model::Model)
    # Return vector of AbstractElement
    vcat(model.frame_elements, model.shell_elements)
end

"""Total number of elements."""
n_elements(model::Model) = model.nFrameElements + model.nShellElements

# For FrameModel/ShellModel compatibility
has_frame_elements(model::FrameModel) = true
has_shell_elements(model::FrameModel) = false
has_frame_elements(model::ShellModel) = false
has_shell_elements(model::ShellModel) = true

# =============================================================================
# make_ids! for unified Model
# =============================================================================

function make_ids!(model::Model)
    make_ids!(model.nodes)
    
    # ID frame elements
    for (i, elem) in enumerate(model.frame_elements)
        elem.elementID = i
    end
    
    # ID shell elements (continue from frame element count)
    offset = model.nFrameElements
    for (i, elem) in enumerate(model.shell_elements)
        elem.elementID = offset + i
    end
    
    make_ids!(model.loads)
end

# Generic make_ids! for all ElementModel types (FrameModel, ShellModel, TrussModel)
function make_ids!(model::ElementModel)
    make_ids!(model.nodes)
    make_ids!(model.elements)
    make_ids!(model.loads)
end


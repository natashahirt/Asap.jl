abstract type AbstractModel end
abstract type ElementModel <: AbstractModel end  # Models with single .elements field

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

# =============================================================================
# FrameModel - For frame elements only (beams, columns)
# =============================================================================

"""
    FrameModel(nodes, elements, loads)

Create a structural model with only frame elements (beams/columns).
This is the original ASAP Model, renamed for clarity.

For mixed frame+shell models, use `Model` instead.
"""
mutable struct FrameModel{E,L} <: ElementModel
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
    Pf::Vector{Float64}
    u::Vector{Float64}
    reactions::Vector{Float64}
    compliance::Float64
    tol::Float64
    processed::Bool
    
    function FrameModel(nodes::Vector{Node}, elements::Vector{E}, loads::Vector{L}) where {E<:FrameElement, L<:AbstractLoad}
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
            zeros(6nnodes),
            0.0,
            1e-6,
            false
        )

        return structure
    end
end

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
            false
        )

        return structure
    end
end

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
    
    # Full constructor: frame + shell elements
    function Model(
        nodes::Vector{Node},
        frame_elements::Vector{E1},
        shell_elements::Vector{E2},
        loads::Vector{L}
    ) where {E1<:FrameElement, E2<:ShellElement, L<:AbstractLoad}
        
        nnodes = length(nodes)
        n_frame = length(frame_elements)
        n_shell = length(shell_elements)

        new(
            nodes,
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
            false
        )
    end
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

# Backward compatibility: allow .elements and .nElements on unified Model
# This maps to frame_elements for code that only uses frame elements
# Provide .elements as alias for .frame_elements for backward compatibility
# Also provide .shells as alias for .shell_elements for unified API
function Base.getproperty(model::Model, name::Symbol)
    if name === :elements
        return getfield(model, :frame_elements)
    elseif name === :nElements
        return getfield(model, :nFrameElements)
    elseif name === :shells
        return getfield(model, :shell_elements)
    elseif name === :nShells
        return getfield(model, :nShellElements)
    else
        return getfield(model, name)
    end
end

# Provide .shells as alias for .elements in ShellModel for unified API
function Base.getproperty(model::ShellModel, name::Symbol)
    if name === :shells
        return getfield(model, :elements)
    else
        return getfield(model, name)
    end
end

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
            false
        )

        return structure
    end
end


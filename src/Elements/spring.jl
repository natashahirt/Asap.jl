#=
Spring Element for Asap.jl
==========================

Grounded spring elements that attach to nodes and provide elastic support.
Useful for modeling:
- Elastic foundations (soil springs)
- Partial fixity conditions
- Spring supports / isolators

Based on FinEtools grounded spring formulations.
=#

using LinearAlgebra: diagm

"""
    Spring(node, k, id=:spring)
    Spring(node; kx=0, ky=0, kz=0, krx=0, kry=0, krz=0, id=:spring)

Create a grounded spring element attached to a node.

Grounded springs connect a node to "ground" with specified stiffness values,
adding diagonal terms to the global stiffness matrix at that node's DOFs.

# Arguments
- `node::Node`: Node to attach spring to
- `k`: Stiffness - either a single value (applied to all translational DOFs),
       a 3-element vector [kx, ky, kz], or a 6-element vector [kx, ky, kz, krx, kry, krz]
- `id::Symbol`: Optional identifier (default: `:spring`)

# Keyword Constructor
- `kx, ky, kz`: Translational stiffnesses (force/length units, e.g., N/m or kip/ft)
- `krx, kry, krz`: Rotational stiffnesses (moment/angle units, e.g., N*m/rad or kip*ft/rad)

# Examples
```julia
# Uniform translational spring (pinned on elastic support)
spring = Spring(node, 1000u"kN/m")

# Different stiffness per translational DOF
spring = Spring(node, [100u"kN/m", 100u"kN/m", 500u"kN/m"])

# Full 6-DOF spring with rotational stiffness
spring = Spring(node; kx=100u"kN/m", ky=100u"kN/m", kz=500u"kN/m", 
                      krz=1000u"kN*m/rad")

# Keyword form (most explicit)
spring = Spring(node; kx=100u"kN/m", ky=100u"kN/m", kz=500u"kN/m")
```

# Notes
- Stiffness values are converted to SI (N/m for translation, N*m/rad for rotation)
- Zero stiffness means that DOF is free (no spring contribution)
- Combine with node fixity: a fixed DOF with a spring is still fixed
- Springs add to the stiffness matrix; they don't replace boundary conditions
"""
mutable struct Spring <: AbstractElement
    node::Node                          # Node to attach spring to
    stiffness::Vector{Float64}          # [kx, ky, kz, krx, kry, krz] in SI units
    elementID::Int64                    # Element index
    globalID::Vector{Int64}             # Element global DOFs (6 for frame node)
    K::Matrix{Float64}                  # 6×6 diagonal stiffness matrix
    id::Symbol                          # Identifier

    # Internal constructor
    function Spring(node::Node, stiff::Vector{Float64}, id::Symbol)
        @assert length(stiff) == 6 "Stiffness vector must have 6 elements"
        K = diagm(stiff)
        new(node, stiff, 0, Vector{Int64}(undef, 6), K, id)
    end
end

# =============================================================================
# Constructors with Unitful quantities
# =============================================================================

# Single scalar: apply to all translational DOFs
function Spring(node::Node, k::Unitful.Quantity{<:Real}, id::Symbol=:spring)
    # Handle both force/length (N/m) and force/length dimension
    k_si = ustrip(Float64, u"N/m", k)
    stiff = [k_si, k_si, k_si, 0.0, 0.0, 0.0]
    Spring(node, stiff, id)
end

# 3-element vector: translational stiffnesses only
function Spring(node::Node, k::AbstractVector{<:Unitful.Quantity}, id::Symbol=:spring)
    if length(k) == 3
        k_si = [ustrip(Float64, u"N/m", ki) for ki in k]
        stiff = [k_si..., 0.0, 0.0, 0.0]
    elseif length(k) == 6
        # First 3 are translational (N/m), last 3 are rotational (N*m/rad)
        k_trans = [ustrip(Float64, u"N/m", k[i]) for i in 1:3]
        k_rot = [ustrip(Float64, u"N*m/rad", k[i]) for i in 4:6]
        stiff = [k_trans..., k_rot...]
    else
        error("Stiffness vector must be length 3 or 6, got $(length(k))")
    end
    Spring(node, stiff, id)
end

# Keyword constructor for explicit per-DOF specification
function Spring(node::Node;
                kx::Unitful.Quantity=0.0u"N/m",
                ky::Unitful.Quantity=0.0u"N/m",
                kz::Unitful.Quantity=0.0u"N/m",
                krx::Unitful.Quantity=0.0u"N*m/rad",
                kry::Unitful.Quantity=0.0u"N*m/rad",
                krz::Unitful.Quantity=0.0u"N*m/rad",
                id::Symbol=:spring)
    stiff = [
        ustrip(Float64, u"N/m", kx),
        ustrip(Float64, u"N/m", ky),
        ustrip(Float64, u"N/m", kz),
        ustrip(Float64, u"N*m/rad", krx),
        ustrip(Float64, u"N*m/rad", kry),
        ustrip(Float64, u"N*m/rad", krz)
    ]
    Spring(node, stiff, id)
end

# =============================================================================
# Dimensionless constructors (assumes SI units)
# =============================================================================

function Spring(node::Node, k::Real, id::Symbol=:spring)
    k_f = Float64(k)
    stiff = [k_f, k_f, k_f, 0.0, 0.0, 0.0]
    Spring(node, stiff, id)
end

function Spring(node::Node, k::AbstractVector{<:Real}, id::Symbol=:spring)
    if length(k) == 3
        stiff = Float64[k..., 0.0, 0.0, 0.0]
    elseif length(k) == 6
        stiff = Float64.(k)
    else
        error("Stiffness vector must be length 3 or 6, got $(length(k))")
    end
    Spring(node, stiff, id)
end

# =============================================================================
# Accessors
# =============================================================================

"""
    stiffness(spring::Spring) -> Vector{Float64}

Get the 6-DOF stiffness vector [kx, ky, kz, krx, kry, krz] in SI units.
"""
spring_stiffness(s::Spring) = s.stiffness

"""
    translational_stiffness(spring::Spring) -> Vector{Float64}

Get translational stiffnesses [kx, ky, kz] in N/m.
"""
translational_stiffness(s::Spring) = s.stiffness[1:3]

"""
    rotational_stiffness(spring::Spring) -> Vector{Float64}

Get rotational stiffnesses [krx, kry, krz] in N*m/rad.
"""
rotational_stiffness(s::Spring) = s.stiffness[4:6]

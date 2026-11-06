# ============================================================================
# SupportLine - internal representation of a support line (forward declaration)
# ============================================================================
# Defined here because _collect_refinement_targets (below) references the type.
# Helper constructors and converters follow later in the file.

"""
Internal struct representing a support line for meshing.
"""
struct _SupportLine
    start_pt::Tuple{Float64, Float64}  # (x, y) in meters
    end_pt::Tuple{Float64, Float64}
end

"""
    _collect_refinement_targets(interior_nodes, support_lines, patches)
        -> Vector{Tuple{Float64,Float64}}

Auto-detect mesh refinement targets from interior features only.
Returns point locations (in meters) where the mesh should be refined.

Targets are:
1. Interior nodes (column connection points)
2. Support-line endpoints (beam connections)
3. Patch centroids (column footprints — stress concentrations)

Boundary conditions (`edge_support_type`) are **not** considered here.
Mesh refinement and DOF fixity are orthogonal concerns.
"""
function _collect_refinement_targets(
    interior_nodes::Vector{Node},
    support_lines::Vector{_SupportLine},
    patches::Vector  # Vector{ShellPatch} — untyped to avoid forward-declaration order
)
    targets = Tuple{Float64, Float64}[]
    
    # 1. Interior nodes (column connection points)
    for node in interior_nodes
        push!(targets, (ustrip(u"m", node.position[1]), 
                        ustrip(u"m", node.position[2])))
    end
    
    # 2. Support line endpoints (beam connections)
    for line in support_lines
        push!(targets, line.start_pt)
        push!(targets, line.end_pt)
    end
    
    # 3. Patch centroids (column footprints — stress concentrations)
    for patch in patches
        push!(targets, patch.center)
    end
    
    # Deduplicate
    return unique(targets)
end

# ============================================================================
# ShellSection - combines thickness and material properties
# ============================================================================

"""
    ShellSection(thickness, E, ν; ρ=0.0u"kg/m^3", name=:section)
    ShellSection(thickness, material::ShellMaterial)

A shell section combining thickness with material properties.

# Examples
```julia
section = ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3")
section = ShellSection(0.15u"m", concrete_material)
```
"""
struct ShellSection
    thickness::Float64  # meters
    E::Float64          # Pa
    ν::Float64
    ρ::Float64          # kg/m³
    name::Symbol
    
    function ShellSection(
        thickness::Quantity,
        E::Quantity,
        ν::Real;
        ρ::Quantity = 0.0u"kg/m^3",
        name::Symbol = :section
    )
        @assert 0.0 < ν < 0.5 "Poisson's ratio must be in (0, 0.5)"
        t = ustrip(u"m", thickness)
        E_val = ustrip(u"Pa", E)
        ρ_val = ustrip(u"kg/m^3", ρ)
        return new(t, E_val, Float64(ν), ρ_val, name)
    end
    
    function ShellSection(thickness::Quantity, material::ShellMaterial)
        t = ustrip(u"m", thickness)
        E_val = ustrip(material.E)
        ρ_val = ustrip(material.ρ)
        return new(t, E_val, material.ν, ρ_val, material.name)
    end
end

# Note: ShellSection presets removed. Create sections inline:
# ShellSection(0.15u"m", 30u"GPa", 0.2; ρ=2400u"kg/m^3", name=:concrete)


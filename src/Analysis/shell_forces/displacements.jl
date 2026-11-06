"""
    ShellDisplacements(elem::ShellElement, u_global::Vector{Float64})

Compute shell displacement at element centroid.
Works for all shell element types (ShellTri3, CompositeShellTri3).
"""
function ShellDisplacements(elem::ShellElement, u_global::Vector{Float64})
    u_local = _extract_local_dofs(elem, u_global)
    
    # Centroid at ξ = η = 1/3: N1 = N2 = N3 = 1/3
    N = 1.0 / 3.0
    
    u = N * (u_local[1] + u_local[7] + u_local[13])
    v = N * (u_local[2] + u_local[8] + u_local[14])
    w = N * (u_local[3] + u_local[9] + u_local[15])
    θx = N * (u_local[4] + u_local[10] + u_local[16])
    θy = N * (u_local[5] + u_local[11] + u_local[17])
    θz = N * (u_local[6] + u_local[12] + u_local[18])
    
    ShellDisplacements(elem, u, v, w, θx, θy, θz)
end

"""Compute displacements for all shell elements in a model."""
function ShellDisplacements(model::ShellModel)
    [ShellDisplacements(elem, model.u) for elem in model.elements]
end

function ShellDisplacements(model::Model)
    [ShellDisplacements(elem, model.u) for elem in model.shell_elements]
end

"""
    max_deflection(disps::Vector{ShellDisplacements})

Find the maximum out-of-plane deflection across all shells.

# Returns
Named tuple: `(max_w, element)` with max deflection value and corresponding element.
"""
function max_deflection(disps::Vector{ShellDisplacements})
    isempty(disps) && return (max_w=0.0, element=nothing)
    max_idx = argmax(abs.(getfield.(disps, :w)))
    (max_w=disps[max_idx].w, element=disps[max_idx].element)
end


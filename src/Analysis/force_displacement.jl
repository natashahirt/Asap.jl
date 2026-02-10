# =============================================================================
# ElementForceAndDisplacement — combined forces + displacements in one pass
# =============================================================================
# Split from force_analysis.jl because this struct depends on both
# ElementInternalForces (force_analysis.jl) and ElementDisplacements
# (displacements.jl), so it must be included after both.

"""
    ElementForceAndDisplacement

Combined internal forces and displacements for a frame element, computed in a
single pass over loads.  Avoids duplicating L / xinc / load iteration that
occurs when `ElementInternalForces` and `ElementDisplacements` are called
separately.

# Fields
- `forces::ElementInternalForces`
- `displacements::ElementDisplacements`
"""
struct ElementForceAndDisplacement
    forces::ElementInternalForces
    displacements::ElementDisplacements
end

"""
    ElementForceAndDisplacement(element, loads; resolution=20)

Compute both internal forces and displacements for `element` from pre-computed
`loads` in a single pass.  Shares `L`, `xinc`, `ElementProps`, and the load
iteration loop.
"""
function ElementForceAndDisplacement(element::AbstractElement, loads::AbstractVector{<:AbstractLoad}; resolution::Int = 20)
    dofs = etype2DOF[typeof(element)]
    ep = ElementProps(element)
    L = ep.L

    rng = range(0, L, resolution)

    # ── Shared: build global displacement vector once (avoids 2 allocs + concat) ──
    _ug = Vector{Float64}(undef, 12)
    _fill_uglobal!(_ug, element)

    # ── Forces: compute Flocal = R * K * uglobal .* dofs ──
    Ku = Vector{Float64}(undef, 12)
    mul!(Ku, element.K, _ug)
    Flocal = Vector{Float64}(undef, 12)
    mul!(Flocal, element.R, Ku)
    @inbounds for i in 1:12
        Flocal[i] *= dofs[i]
    end

    Pstart  = -Flocal[1]
    Vystart =  Flocal[2]
    Mystart =  Flocal[6]
    Vzstart =  Flocal[3]
    Mzstart = -Flocal[5]

    P  = fill(Pstart, resolution)
    My = Vector{Float64}(undef, resolution)
    Vy = fill(Vystart, resolution)
    Mz = Vector{Float64}(undef, resolution)
    Vz = fill(Vzstart, resolution)
    @inbounds for i in 1:resolution
        My[i] = Vystart * rng[i] - Mystart
        Mz[i] = Vzstart * rng[i] - Mzstart
    end

    # ── Displacements: compute ulocal = R * uglobal .* dofs (reuses _ug) ──
    _ul = Vector{Float64}(undef, 12)
    mul!(_ul, element.R, _ug)
    @inbounds for i in 1:12
        _ul[i] *= dofs[i]
    end

    uX1 = _ul[1];  uX2 = _ul[7]
    uY1 = _ul[2];  uY2 = _ul[6];  uY3 = _ul[8];  uY4 = _ul[12]
    uZ1 = _ul[3];  uZ2 = -_ul[5]; uZ3 = _ul[9];  uZ4 = -_ul[11]

    Dx = Vector{Float64}(undef, resolution)
    Dy = Vector{Float64}(undef, resolution)
    Dz = Vector{Float64}(undef, resolution)

    @inbounds for i in 1:resolution
        x = rng[i]
        xL = x / L

        # Axial: linear interpolation
        Dx[i] = (1 - xL) * uX1 + xL * uX2

        # Transverse (Hermite): N(x,L) * u
        n1 = 1 - 3xL^2 + 2xL^3
        n2 = x * (1 - xL)^2
        n3 = 3xL^2 - 2xL^3
        n4 = x^2/L * (-1 + xL)

        Dy[i] = n1*uY1 + n2*uY2 + n3*uY3 + n4*uY4
        Dz[i] = n1*uZ1 + n2*uZ2 + n3*uZ3 + n4*uZ4
    end

    # ── Single pass over loads: accumulate both forces and displacements ──
    for load in loads
        accumulate_force!(load, rng, P, My, Vy, Mz, Vz, ep)
        accumulatedisp!(load, rng, Dy, Dz, ep)
    end

    xinc = collect(rng)
    forces = ElementInternalForces(element, resolution, xinc, P, My, Vy, Mz, Vz)

    # ── Build displacement matrices ──
    D = Matrix{Float64}(undef, 3, resolution)
    @inbounds for j in 1:resolution
        D[1,j] = Dx[j]; D[2,j] = Dy[j]; D[3,j] = Dz[j]
    end

    lx, ly, lz = element.LCS[1], element.LCS[2], element.LCS[3]
    Dglobal = Matrix{Float64}(undef, 3, resolution)
    @inbounds for j in 1:resolution
        dx, dy, dz = Dx[j], Dy[j], Dz[j]
        Dglobal[1,j] = dx*lx[1] + dy*ly[1] + dz*lz[1]
        Dglobal[2,j] = dx*lx[2] + dy*ly[2] + dz*lz[2]
        Dglobal[3,j] = dx*lx[3] + dy*ly[3] + dz*lz[3]
    end

    px = ustrip(u"m", element.nodeStart.position[1])
    py = ustrip(u"m", element.nodeStart.position[2])
    pz = ustrip(u"m", element.nodeStart.position[3])
    ax_x, ax_y, ax_z = lx[1], lx[2], lx[3]
    basepoints = Matrix{Float64}(undef, 3, resolution)
    @inbounds for j in 1:resolution
        basepoints[1,j] = px + ax_x * rng[j]
        basepoints[2,j] = py + ax_y * rng[j]
        basepoints[3,j] = pz + ax_z * rng[j]
    end

    disp = ElementDisplacements(element, resolution, xinc, D, Dglobal, basepoints)

    return ElementForceAndDisplacement(forces, disp)
end

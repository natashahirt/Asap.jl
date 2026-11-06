# ============================================================================
# Principal Values and Design Quantities
# ============================================================================

"""
    principal_moments(sif::ShellInternalForces)

Compute principal bending moments and their directions.

# Returns
- `(M1=..., M2=..., θ=...)`: NamedTuple with principal moments and angle (radians)
"""
function principal_moments(sif::ShellInternalForces)
    Mx, My, Mxy = sif.Mxx, sif.Myy, sif.Mxy
    
    M_avg = (Mx + My) / 2
    R = sqrt(((Mx - My) / 2)^2 + Mxy^2)
    
    M1 = M_avg + R
    M2 = M_avg - R
    θ = 0.5 * atan(2 * Mxy, Mx - My)
    
    return (M1=M1, M2=M2, θ=θ)
end

"""
    principal_forces(sif::ShellInternalForces)

Compute principal membrane forces and their directions.

# Returns
- `(N1=..., N2=..., θ=...)`: NamedTuple with principal forces and angle (radians)
"""
function principal_forces(sif::ShellInternalForces)
    Nx, Ny, Nxy = sif.Nxx, sif.Nyy, sif.Nxy
    
    N_avg = (Nx + Ny) / 2
    R = sqrt(((Nx - Ny) / 2)^2 + Nxy^2)
    
    N1 = N_avg + R
    N2 = N_avg - R
    θ = 0.5 * atan(2 * Nxy, Nx - Ny)
    
    return (N1=N1, N2=N2, θ=θ)
end

"""
    von_mises_stress(sif::ShellInternalForces, z::Float64, t::Float64)

Compute von Mises stress at a given z-coordinate through the thickness.

# Arguments
- `z`: Distance from mid-plane (positive = top face)
- `t`: Total thickness

# Returns
von Mises stress in Pa
"""
function von_mises_stress(sif::ShellInternalForces, z::Float64, t::Float64)
    σxx = sif.Nxx / t + 12 * sif.Mxx * z / t^3
    σyy = sif.Nyy / t + 12 * sif.Myy * z / t^3
    τxy = sif.Nxy / t + 12 * sif.Mxy * z / t^3
    
    sqrt(σxx^2 - σxx*σyy + σyy^2 + 3*τxy^2)
end

"""
    max_surface_stresses(sif::ShellInternalForces, t::Float64)

Compute maximum stresses on top (+t/2) and bottom (-t/2) surfaces.

# Returns
Named tuple: `(top=(σxx, σyy, τxy), bottom=(σxx, σyy, τxy))`
"""
function max_surface_stresses(sif::ShellInternalForces, t::Float64)
    σxx_top = sif.Nxx / t + 6 * sif.Mxx / t^2
    σyy_top = sif.Nyy / t + 6 * sif.Myy / t^2
    τxy_top = sif.Nxy / t + 6 * sif.Mxy / t^2
    
    σxx_bot = sif.Nxx / t - 6 * sif.Mxx / t^2
    σyy_bot = sif.Nyy / t - 6 * sif.Myy / t^2
    τxy_bot = sif.Nxy / t - 6 * sif.Mxy / t^2
    
    (top = (σxx=σxx_top, σyy=σyy_top, τxy=τxy_top),
     bottom = (σxx=σxx_bot, σyy=σyy_bot, τxy=τxy_bot))
end


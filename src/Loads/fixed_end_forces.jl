"""
    q_local(load::PointLoad)

Equivalent fixed end forces for a point load.
"""
function q_local(load::PointLoad)
    #values
    LCS = load.element.LCS
    l = to_meters(load.element.length)

    # load vector in LCS (strip units from load.value)
    value_stripped = [to_newtons(v) for v in load.value]
    plocal = load.element.R[1:3, 1:3] * value_stripped .* LCS

    #axial end forces
    ax1 = ax2 = - dot(plocal[1], LCS[1]) / 2

    #position of load from start node
    a = load.position * l
    b = l - a #remainder

    # perpendicular load in local Y
    py = - dot(plocal[2], LCS[2])

    # moments in local Z
    mz1 = py * b^2 * a / l^2 
    # mz2 = -py * a^2 * b / l^2
    mz2 = -py * a^2 * b / l^2

    # shear in local Y
    vy1 = py * b^2 / l^3 * (3a + b)
    vy2 = py * a^2 / l^3 * (a + 3b)

    # perpendicular load in local Z
    pz = - dot(plocal[3], LCS[3])

    # moments in local Y
    my1 = -pz * b^2 * a / l^2
    my2 = pz * a^2 * b / l^2

    # shear in local Z
    vz1 = pz * b^2 / l^3 * (3a + b)
    vz2 = pz * a^2 / l^3 * (a + 3b)

    return [ax1, vy1, vz1, 0., my1, mz1, ax2, vy2, vz2, 0., my2, mz2]
end

"""
    q_local(load::LineLoad)

Equivalent fixed end forces for a line load
"""
function q_local(load::LineLoad)
    LCS = load.element.LCS
    l = to_meters(load.element.length)

    # load vector in LCS (strip units from load.value - it's force/length N/m)
    value_stripped = [to_newtons_per_meter(v) for v in load.value]
    plocal = load.element.R[1:3, 1:3] * value_stripped .* LCS

    #axial end forces
    ax1 = ax2 = - dot(plocal[1], LCS[1]) * l / 2

    #perpendicular load in local y
    py = - dot(plocal[2], LCS[2])

    vy1 = vy2 = py * l / 2 #shears in Y
    mz1 = py * l^2 / 12 #moment 1 in Z
    mz2 = -mz1 #moment 2 in Z


    # perpendicular load in local z
    pz = -dot(plocal[3], LCS[3])

    vz1 = vz2 = pz * l / 2 #shears
    my1 = -pz * l^2 / 12 #moment 1 in Y
    my2 = -my1 #moment 2 in y

    return [ax1, vy1, vz1, 0., my1, mz1, ax2, vy2, vz2, 0., my2, mz2]
end

"""
    q_local(load::LineLoad)

Equivalent fixed end forces for a gravity (line) load
"""
function q_local(load::GravityLoad)

    LCS = load.element.LCS
    l = to_meters(load.element.length)
    value = [0., 0., -1.] .* to_kg_per_m3(load.element.section.ρ) .* to_meters_squared(load.element.section.A) .* to_m_per_s2(load.factor)

    # load vector in LCS
    plocal = load.element.R[1:3, 1:3] * value .* LCS

    #axial end forces
    ax1 = ax2 = - dot(plocal[1], LCS[1]) * l / 2

    #perpendicular load in local y
    py = - dot(plocal[2], LCS[2])

    vy1 = vy2 = py * l / 2 #shears in Y
    mz1 = py * l^2 / 12 #moment 1 in Z
    mz2 = -mz1 #moment 2 in Z

    # perpendicular load in local z
    pz = -dot(plocal[3], LCS[3])

    vz1 = vz2 = pz * l / 2 #shears
    my1 = -pz * l^2 / 12 #moment 1 in Y
    my2 = -my1 #moment 2 in y

    return [ax1, vy1, vz1, 0., my1, mz1, ax2, vy2, vz2, 0., my2, mz2]
end

function q(load::ElementLoad{FixedFixed})

    return q_local(load)

end

function q(load::ElementLoad{FixedFree})
    #length of element
    factor = 3 / 2 / to_meters(load.element.length)
    #fixed end components
    FAb, FSby, FSbz, FTb, FMby, FMbz, FAe, FSey, FSez, FTe, FMey, FMez = q_local(load)

    #modified fixed end forces
    return [FAb, 
        FSby - factor*FMez,
        FSbz + factor*FMey,
        FTb + FTe,
        FMby - 1/2 * FMey,
        FMbz - 1/2 * FMez,
        FAe,
        FSey + factor*FMez,
        FSez - factor*FMey,
        0,
        0,
        0]
end

function q(load::ElementLoad{FreeFixed})
    #length of element
    factor = 3 / 2 / to_meters(load.element.length)
    #fixed end components
    FAb, FSby, FSbz, FTb, FMby, FMbz, FAe, FSey, FSez, FTe, FMey, FMez = q_local(load)

    #modified fixed end forces
    return [FAb, #axial beginning
        FSby - factor*FMbz, #shear local Y
        FSbz + factor*FMby, #shear local Z
        0, #torsion beginning
        0, #moment local Y
        0, #moment local Z
        FAe, # axial end
        FSey + factor*FMbz, #shear local Y
        FSez - factor*FMby, #shear local Z
        FTb + FTe, #torsion end
        FMey - 1/2*FMby, #moment local Y
        FMez - 1/2*FMbz] #moment local Z
end

function q(load::ElementLoad{FreeFree})
    #length of element
    factor = 1 / to_meters(load.element.length)
    #fixed end components
    FAb, FSby, FSbz, FTb, FMby, FMbz, FAe, FSey, FSez, FTe, FMey, FMez = q_local(load)

    #modified fixed end forces
    return [FAb, 
        FSby - factor*(FMbz + FMez),
        FSbz + factor*(FMby + FMey),
        FTb,
        0,
        0,
        FAe,
        FSey + factor*(FMbz + FMez),
        FSez - factor*(FMby + FMey),
        FTe,
        0,
        0]
end

function q(load::ElementLoad{Joist})
    #length of element
    factor = 1 / to_meters(load.element.length)
    #fixed end components
    FAb, FSby, FSbz, FTb, FMby, FMbz, FAe, FSey, FSez, FTe, FMey, FMez = q_local(load)

    #modified fixed end forces
    return [FAb, 
        FSby - factor*(FMbz + FMez),
        FSbz + factor*(FMby + FMey),
        FTb,
        0,
        0,
        FAe,
        FSey + factor*(FMbz + FMez),
        FSez - factor*(FMby + FMey),
        FTe,
        0,
        0]
end

# =============================================================================
# TributaryLoad - Piecewise Linear Distributed Load
# =============================================================================

"""
    q_local(load::TributaryLoad)

Equivalent fixed end forces for a piecewise-linear tributary load.

Uses exact analytical integration for each linear segment.
"""
function q_local(load::TributaryLoad)
    LCS = load.element.LCS
    L = to_meters(load.element.length)
    
    # Get intensities (N/m) at each breakpoint
    w_vals = intensities(load)
    
    # Build global load direction vector and transform to local
    dir_global = collect(load.direction)
    dir_local = load.element.R[1:3, 1:3] * dir_global .* LCS
    
    # Initialize FEM accumulator [ax1, vy1, vz1, 0, my1, mz1, ax2, vy2, vz2, 0, my2, mz2]
    fem = zeros(12)
    
    # Sum contributions from each linear segment
    for i in 1:(length(load.positions) - 1)
        s_a = load.positions[i]
        s_b = load.positions[i + 1]
        w_a = w_vals[i]
        w_b = w_vals[i + 1]
        
        # Skip zero-length or zero-load segments
        (s_b - s_a) < 1e-12 && continue
        (abs(w_a) + abs(w_b)) < 1e-12 && continue
        
        # Compute FEM contribution from this linear segment
        fem .+= _fem_linear_segment(L, s_a, s_b, w_a, w_b, dir_local, LCS)
    end
    
    return fem
end

"""
Compute fixed-end forces for a linear load segment.

Arguments:
- L: beam length (m)
- s_a, s_b: normalized positions [0,1] of segment start/end
- w_a, w_b: load intensities (N/m) at start/end
- dir_local: load direction in local coordinates
- LCS: local coordinate system vectors
"""
function _fem_linear_segment(L::Float64, s_a::Float64, s_b::Float64, 
                             w_a::Float64, w_b::Float64,
                             dir_local::Vector{Vector{Float64}}, 
                             LCS::Vector{Vector{Float64}})
    # Absolute positions
    a = s_a * L
    b = s_b * L
    ℓ = b - a  # segment length
    
    # Decompose into uniform + triangular:
    # w(x) = w_a + (w_b - w_a)*(x - a)/ℓ
    #      = w_a + Δw*(x - a)/ℓ  where Δw = w_b - w_a
    w_uniform = w_a
    Δw = w_b - w_a
    
    # Fixed-end forces for uniform load w over [a, b]
    # Using exact integrals of point-load FEM formulas
    R_A_u, R_B_u, M_A_u, M_B_u = _fem_uniform_partial(L, a, b, w_uniform)
    
    # Fixed-end forces for triangular load (0 at a, Δw at b)
    R_A_t, R_B_t, M_A_t, M_B_t = _fem_triangular_partial(L, a, b, Δw)
    
    # Total reactions and moments
    R_A = R_A_u + R_A_t
    R_B = R_B_u + R_B_t
    M_A = M_A_u + M_A_t
    M_B = M_B_u + M_B_t
    
    # Distribute to FEM vector based on load direction
    fem = zeros(12)
    
    # Axial component (in local X)
    px = -dot(dir_local[1], LCS[1])
    fem[1] += px * (R_A + R_B) / 2  # ax1 (simplified: half total)
    fem[7] += px * (R_A + R_B) / 2  # ax2
    
    # Transverse in local Y
    py = -dot(dir_local[2], LCS[2])
    fem[2] += py * R_A         # vy1 (shear at start)
    fem[8] += py * R_B         # vy2 (shear at end)
    fem[6] += py * M_A         # mz1 (moment about Z at start)
    fem[12] += py * M_B        # mz2 (moment about Z at end)
    
    # Transverse in local Z
    pz = -dot(dir_local[3], LCS[3])
    fem[3] += pz * R_A         # vz1
    fem[9] += pz * R_B         # vz2
    fem[5] += -pz * M_A        # my1 (note sign convention)
    fem[11] += -pz * M_B       # my2
    
    return fem
end

"""
Fixed-end forces for uniform load w over partial span [a, b] on beam of length L.
Returns (R_A=..., R_B=..., M_A=..., M_B=...) where positive M is counterclockwise.
"""
function _fem_uniform_partial(L::Float64, a::Float64, b::Float64, w::Float64)
    abs(w) < 1e-15 && return (R_A=0.0, R_B=0.0, M_A=0.0, M_B=0.0)
    
    # Integrals for fixed-fixed beam FEMs
    # R_A = ∫[a,b] w*(L-x)²*(L+2x)/L³ dx
    # R_B = ∫[a,b] w*x²*(3L-2x)/L³ dx
    # M_A = -∫[a,b] w*x*(L-x)²/L² dx
    # M_B = ∫[a,b] w*x²*(L-x)/L² dx
    
    L2 = L^2
    L3 = L^3
    
    # Evaluate antiderivatives at b and a
    # For R_A: F(x) = (L-x)³*(L+2x)/(-3L³) + 2(L-x)⁴/(12L³) ... using integration by parts
    # Let's use direct polynomial integration instead
    
    # Expand (L-x)²*(L+2x) = (L²-2Lx+x²)(L+2x) = L³+2L²x-2L²x-4Lx²+Lx²+2x³ = L³-3Lx²+2x³
    # Wait, let me redo: (L-x)² = L²-2Lx+x², (L+2x) = L+2x
    # (L²-2Lx+x²)(L+2x) = L³ + 2L²x - 2L²x - 4Lx² + Lx² + 2x³ = L³ - 3Lx² + 2x³
    
    # R_A = w/L³ * ∫[a,b] (L³ - 3Lx² + 2x³) dx = w/L³ * [L³x - Lx³ + x⁴/2]|[a,b]
    I_RA(x) = L3*x - L*x^3 + x^4/2
    R_A = w/L3 * (I_RA(b) - I_RA(a))
    
    # For R_B: x²*(3L-2x) = 3Lx² - 2x³
    # R_B = w/L³ * ∫[a,b] (3Lx² - 2x³) dx = w/L³ * [Lx³ - x⁴/2]|[a,b]
    I_RB(x) = L*x^3 - x^4/2
    R_B = w/L3 * (I_RB(b) - I_RB(a))
    
    # For M_A: -x*(L-x)² = -x*(L² - 2Lx + x²) = -L²x + 2Lx² - x³
    # M_A = -w/L² * ∫[a,b] (L²x - 2Lx² + x³) dx = -w/L² * [L²x²/2 - 2Lx³/3 + x⁴/4]|[a,b]
    I_MA(x) = L2*x^2/2 - 2L*x^3/3 + x^4/4
    M_A = -w/L2 * (I_MA(b) - I_MA(a))
    
    # For M_B: x²*(L-x) = Lx² - x³
    # M_B = w/L² * ∫[a,b] (Lx² - x³) dx = w/L² * [Lx³/3 - x⁴/4]|[a,b]
    I_MB(x) = L*x^3/3 - x^4/4
    M_B = w/L2 * (I_MB(b) - I_MB(a))
    
    return (R_A=R_A, R_B=R_B, M_A=M_A, M_B=M_B)
end

"""
Fixed-end forces for triangular load (0 at a, w_peak at b) over [a, b] on beam of length L.
Returns (R_A=..., R_B=..., M_A=..., M_B=...) where positive M is counterclockwise.
"""
function _fem_triangular_partial(L::Float64, a::Float64, b::Float64, w_peak::Float64)
    abs(w_peak) < 1e-15 && return (R_A=0.0, R_B=0.0, M_A=0.0, M_B=0.0)
    
    ℓ = b - a
    ℓ < 1e-12 && return (R_A=0.0, R_B=0.0, M_A=0.0, M_B=0.0)
    
    # Load function: w(x) = w_peak * (x - a) / ℓ  for a ≤ x ≤ b
    # Need to integrate: w_peak/ℓ * (x - a) * [FEM kernel]
    
    L2 = L^2
    L3 = L^3
    
    # For R_A: ∫[a,b] (w_peak/ℓ)*(x-a)*(L-x)²*(L+2x)/L³ dx
    # Expand (x-a)*(L-x)²*(L+2x) = (x-a)*(L³ - 3Lx² + 2x³)
    # = L³x - 3Lx³ + 2x⁴ - aL³ + 3aLx² - 2ax³
    # = -aL³ + L³x + 3aLx² - 3Lx³ - 2ax³ + 2x⁴
    I_RA(x) = -a*L3*x + L3*x^2/2 + a*L*x^3 - 3L*x^4/4 - 2a*x^4/4 + 2x^5/5
    R_A = (w_peak/ℓ)/L3 * (I_RA(b) - I_RA(a))
    
    # For R_B: ∫[a,b] (w_peak/ℓ)*(x-a)*x²*(3L-2x)/L³ dx
    # (x-a)*x²*(3L-2x) = (x-a)*(3Lx² - 2x³) = 3Lx³ - 2x⁴ - 3aLx² + 2ax³
    I_RB(x) = 3L*x^4/4 - 2x^5/5 - a*L*x^3 + a*x^4/2
    R_B = (w_peak/ℓ)/L3 * (I_RB(b) - I_RB(a))
    
    # For M_A: -∫[a,b] (w_peak/ℓ)*(x-a)*x*(L-x)²/L² dx
    # (x-a)*x*(L-x)² = (x-a)*x*(L² - 2Lx + x²) = (x-a)*(L²x - 2Lx² + x³)
    # = L²x² - 2Lx³ + x⁴ - aL²x + 2aLx² - ax³
    I_MA(x) = L2*x^3/3 - 2L*x^4/4 + x^5/5 - a*L2*x^2/2 + 2a*L*x^3/3 - a*x^4/4
    M_A = -(w_peak/ℓ)/L2 * (I_MA(b) - I_MA(a))
    
    # For M_B: ∫[a,b] (w_peak/ℓ)*(x-a)*x²*(L-x)/L² dx
    # (x-a)*x²*(L-x) = (x-a)*(Lx² - x³) = Lx³ - x⁴ - aLx² + ax³
    I_MB(x) = L*x^4/4 - x^5/5 - a*L*x^3/3 + a*x^4/4
    M_B = (w_peak/ℓ)/L2 * (I_MB(b) - I_MB(a))
    
    return (R_A=R_A, R_B=R_B, M_A=M_A, M_B=M_B)
end
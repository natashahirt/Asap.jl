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
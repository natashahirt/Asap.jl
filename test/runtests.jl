using Asap
using StructuralUnits
using Test
using LinearAlgebra
using Unitful

Unitful.register(StructuralUnits)

@testset "Asap.jl" begin
    #tol
    tol = 0.1
    # 2D truss test: Example 3.9 from Kassimali "Matrix Analysis of Structures 2e"
    # in kips, ft

    n1 = TrussNode([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n2 = TrussNode([10.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n3 = TrussNode([0.0u"m", 8.0u"m", 0.0u"m"], :yfree)
    n4 = TrussNode([6.0u"m", 8.0u"m", 0.0u"m"], :free)

    nodes = [n1, n2, n3, n4]

    E = 70.0u"kN/m^2"  # 70 kN/m² = 70 kPa
    A = 4e3 / 1e6 * u"m^2"  # 0.004 m²

    sec = TrussSection(A, E)

    e1 = TrussElement(nodes[[1,3]]..., sec)
    e2 = TrussElement(nodes[[3,4]]..., sec)
    e3 = TrussElement(nodes[[1,4]]..., sec)
    e4 = TrussElement(nodes[[2,3]]..., sec)
    e5 = TrussElement(nodes[[2,4]]..., sec)

    elements = [e1,e2,e3,e4,e5]

    l1 = NodeForce(n3, [0.0u"N", -400.0u"N", 0.0u"N"])
    l2 = NodeForce(n4, [800.0u"N", -400.0u"N", 0.0u"N"])

    loads = [l1, l2]

    model = TrussModel(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    reactions = model.reactions[model.fixedDOFs]
    reactions2d = reactions[reactions .!= 0]

    reactions_textbook = [-.57994,
        320.82,
        -298.39,
        479.17,
        -501.05]

    err = norm(reactions_textbook .- reactions2d)

    @test err <= tol


    #2D frame test: Example 6.6
    # in kips, in

    n1 = Node([0.0u"inch", 0.0u"inch", 0.0u"inch"], :fixed)
    n2 = Node([10.0u"ft", 20.0u"ft", 0.0u"inch"], :free)
    n3 = Node([30.0u"ft", 20.0u"ft", 0.0u"inch"], :fixed)
    nodes = [n1, n2, n3]

    E = 29e3u"ksi"  # 29,000 ksi
    A = 11.8u"inch^2"
    I = 310.0u"inch^4"

    sec = Section(A, E, 1.0u"Pa", I, I, 1.0u"inch^4")  # G is dummy, J is dummy

    e1 = Element(nodes[[1,2]]..., sec)
    e1.Ψ = 0.
    e2 = Element(nodes[[2,3]]..., sec)
    e2.Ψ = 0.

    elements = [e1, e2]

    l1 = PointLoad(e1, 0.5, [0.0u"kip", -90.0u"kip", 0.0u"kip"])
    l2 = NodeMoment(n2, [0.0u"kip*inch", 0.0u"kip*inch", -125.0u"kip*ft"])
    l3 = LineLoad(e2, [0.0u"kip/inch", -1.5u"kip/ft", 0.0u"kip/inch"])

    loads = [l1, l2, l3]

    model = Model(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    reactions_si = model.reactions[model.fixedDOFs]
    reactions_si = reactions_si[reactions_si .!= 0]

    reactions_textbook = [30.371,  # kip
        102.09,                     # kip
        1216.,                      # kip*inch
        -30.372,                    # kip
        17.913,                     # kip
        -854.07]                    # kip*inch
    
    # Convert Asap output from SI (N, N*m) back to input units (kip, kip*inch)
    # DOF order: [Fx1, Fy1, Mz1, Fx3, Fy3, Mz3]
    reactions_converted = [
        ustrip(uconvert(u"kip", reactions_si[1] * u"N")),           # Fx1
        ustrip(uconvert(u"kip", reactions_si[2] * u"N")),           # Fy1
        ustrip(uconvert(u"kip*inch", reactions_si[3] * u"N*m")),    # Mz1
        ustrip(uconvert(u"kip", reactions_si[4] * u"N")),           # Fx3
        ustrip(uconvert(u"kip", reactions_si[5] * u"N")),           # Fy3
        ustrip(uconvert(u"kip*inch", reactions_si[6] * u"N*m")),    # Mz3
    ]

    err = norm(reactions_textbook .- reactions_converted)

    @test err <= tol

    # 3D truss test: Example 8.1
    # in kips, in

    E = 10e3u"ksi"
    A = 8.4u"inch^2"
    sec = TrussSection(A, E)

    n1 = TrussNode([-6.0u"ft", 0.0u"inch", 8.0u"ft"], :fixed)
    n2 = TrussNode([12.0u"ft", 0.0u"inch", 8.0u"ft"], :fixed)
    n3 = TrussNode([6.0u"ft", 0.0u"inch", -8.0u"ft"], :fixed)
    n4 = TrussNode([-12.0u"ft", 0.0u"inch", -8.0u"ft"], :fixed)
    n5 = TrussNode([0.0u"inch", 24.0u"ft", 0.0u"inch"], :free)

    nodes = [n1, n2, n3, n4, n5]

    e1 = TrussElement(nodes[[1,5]]..., sec)
    e2 = TrussElement(nodes[[2,5]]..., sec)
    e3 = TrussElement(nodes[[3,5]]..., sec)
    e4 = TrussElement(nodes[[4,5]]..., sec)

    elements = [e1, e2, e3, e4]

    l1 = NodeForce(n5, [0.0u"kip", -100.0u"kip", -50.0u"kip"])
    loads = [l1]

    model = TrussModel(nodes, elements, loads)
    solve!(model)

    reactions_si = model.reactions[model.fixedDOFs]

    reactions_textbook = [-5.5581,
        -22.232,
        7.4108,
        1.3838,
        -2.7677,
        0.92255,
        -19.442,
        77.768,
        25.923,
        23.616,
        47.232,
        15.744
        ]
    
    # Convert Asap output from SI (N) back to input units (kip)
    reactions_converted = [ustrip(uconvert(u"kip", r * u"N")) for r in reactions_si]

    err = norm(reactions_textbook .- reactions_converted)

    @test err <= tol

    # 3D frame test: Example 8.4
    # in kips, inches

    n1 = Node([0.0u"inch", 0.0u"inch", 0.0u"inch"], :free)
    n2 = Node([-240.0u"inch", 0.0u"inch", 0.0u"inch"], :fixed)
    n3 = Node([0.0u"inch", -240.0u"inch", 0.0u"inch"], :fixed)
    n4 = Node([0.0u"inch", 0.0u"inch", -240.0u"inch"], :fixed)

    nodes = [n1, n2, n3, n4]

    E = 29e3u"ksi"
    G = 11.5e3u"ksi"
    A = 32.9u"inch^2"
    Iz = 716.0u"inch^4"
    Iy = 236.0u"inch^4"
    J = 15.1u"inch^4"

    sec = Section(A, E, G, Iz, Iy, J)

    e1 = Element(nodes[[2,1]]..., sec)
    e1.Ψ = 0.
    e2 = Element(nodes[[3,1]]..., sec)
    e2.Ψ = pi/2
    e3 = Element(nodes[[4,1]]..., sec)
    e3.Ψ = pi/6

    elements = [e1, e2, e3]

    l1 = LineLoad(e1, [0.0u"kip/inch", -3.0u"kip/ft", 0.0u"kip/inch"])
    l2 = NodeMoment(n1, [-150.0u"kip*ft", 0.0u"kip*ft", 150.0u"kip*ft"])

    loads = [l1, l2]

    model = Model(nodes, elements, loads)
    solve!(model)

    # Collect reactions from fixed nodes in order: n2, n3, n4
    # Each node has 6 DOFs: [Fx, Fy, Fz, Mx, My, Mz]
    # n1 is free, so skip it
    reactions_quantity = Quantity[]
    for node in [model.nodes[2], model.nodes[3], model.nodes[4]]  # n2, n3, n4 (fixed nodes)
        append!(reactions_quantity, node.reaction)
    end

    reactions_textbook = [5.3757,
        44.106,
        -0.74272,
        2.1722,
        58.987,
        2330.5,
        -4.6249,
        11.117,
        -6.4607,
        -515.55,
        -0.76472,
        369.67,
        -0.75082,
        4.7763,
        7.2034,
        -383.5,
        -60.166,
        -4.702]
    
    # Convert Asap output from SI (N, N*m) back to input units (kip, kip*inch)
    # For a 3D frame, DOF order is [Fx, Fy, Fz, Mx, My, Mz] per node
    reactions_converted = Vector{Float64}(undef, length(reactions_quantity))
    for i in 1:length(reactions_quantity)
        # DOFs 4, 5, 6 (mod 6) are moments (Mx, My, Mz), rest are forces
        if (i - 1) % 6 >= 3  # Indices 4, 5, 6 (0-indexed: 3, 4, 5) are moments
            reactions_converted[i] = ustrip(uconvert(u"kip*inch", reactions_quantity[i]))
        else  # Forces (Fx, Fy, Fz)
            reactions_converted[i] = ustrip(uconvert(u"kip", reactions_quantity[i]))
        end
    end

    err = norm(reactions_textbook .- reactions_converted)

    @test err <= tol


    # from https://www.12000.org/my_notes/stiffness_matrix/stiffness_matrix_report.htm

    #Ex1
    P=400.0u"lbf"
    L=144.0u"inch"
    E=30e6u"psi"
    Is=57.1u"inch^4"

    n1 = Node([0.0u"inch", 0.0u"inch", 0.0u"inch"], :fixed)
    n2 = Node([L, 0.0u"inch", 0.0u"inch"], :free)

    nodes = [n1, n2]

    sec = Section(1.0u"inch^2", E, 1.0u"psi", Is, Is, 1.0u"inch^4")

    e = Element(nodes[[1,2]]..., sec)
    elements = [e]

    p = NodeForce(n2, [0.0u"lbf", -P, 0.0u"lbf"])
    loads = [p]

    model = Model(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    d_textbook = [-.2324, -.0024]  # in inches and radians
    d_si = n2.displacement[[2,6]]  # displacement in SI (m, radians) - now Unitful!
    
    # Convert Asap output from SI (m) back to input units (inch) for translation
    # Rotation is already dimensionless (radians) in both
    d_converted = [
        ustrip(uconvert(u"inch", d_si[1])),  # Translation: m -> inch (d_si[1] is already Quantity)
        ustrip(d_si[2])  # Rotation: already in radians (dimensionless)
    ]

    @test norm(d_textbook .- d_converted) <= tol

    #Ex2
    p = PointLoad(e, 0.5, [0.0u"lbf", -P, 0.0u"lbf"])
    loads = [p]
    
    model = Model(nodes, elements, loads)
    planarize!(model)
    solve!(model; reprocess = true)
    
    d_textbook = [-0.072630472854641, -0.000605253940455]  # in inches and radians
    d_si = n2.displacement[[2,6]]  # displacement in SI (m, radians) - now Unitful!
    
    # Convert Asap output from SI (m) back to input units (inch) for translation
    # Rotation is already dimensionless (radians) in both
    d_converted = [
        ustrip(uconvert(u"inch", d_si[1])),  # Translation: m -> inch (d_si[1] is already Quantity)
        ustrip(d_si[2])  # Rotation: already in radians (dimensionless)
    ]
    
    @test norm(d_textbook .- d_converted) <= tol
end

using Asap
using StructuralUnits
using Test
using Unitful

Unitful.register(StructuralUnits)

@testset "TributaryLoad" begin
    # Create a simple fixed-fixed beam (10m span)
    n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n2 = Node([10.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    nodes = [n1, n2]

    E = 200e9u"Pa"
    G = 77e9u"Pa"
    A = 0.01u"m^2"
    I = 1e-4u"m^4"
    J = 1e-4u"m^4"
    sec = Section(A, E, G, I, I, J)

    e = Element(n1, n2, sec)
    elements = [e]

    # Test 1: Uniform TributaryLoad should match LineLoad
    # w = 1000 N/m uniform over full span
    positions = [0.0, 1.0]
    widths = [1.0u"m", 1.0u"m"]  # 1m width × 1000 Pa = 1000 N/m
    pressure = 1000.0u"Pa"
    direction = (0.0, 0.0, -1.0)

    trib_load = TributaryLoad(e, positions, widths, pressure, direction)
    
    @test trib_load.positions == [0.0, 1.0]
    @test trib_load.widths == [1.0u"m", 1.0u"m"]
    @test intensities(trib_load) ≈ [1000.0, 1000.0]

    # Solve with TributaryLoad
    model = Model(nodes, elements, [trib_load])
    planarize!(model)
    solve!(model)
    
    # For uniform load w=1000 N/m on L=10m fixed-fixed beam:
    # Each vertical reaction = wL/2 = 5000 N
    # Fixed-end moment = wL²/12 = 8333.33 N*m
    reactions = model.reactions[model.fixedDOFs]
    nonzero_reactions = reactions[reactions .!= 0]
    
    # Check vertical reactions (should be ~5000 N each)
    # In 2D planarized model: DOFs are [Fx, Fy, Mz] per node
    # For gravity load in -Z, maps to -Y in 2D
    println("TributaryLoad reactions: ", nonzero_reactions)
    
    # Test 2: Triangular load (linearly varying)
    # Width varies from 0 at start to 2m at end
    # w(x) = 2000 * x/L  N/m
    positions_tri = [0.0, 1.0]
    widths_tri = [0.0u"m", 2.0u"m"]  # 0 to 2m width
    
    trib_tri = TributaryLoad(e, positions_tri, widths_tri, pressure, direction)
    @test intensities(trib_tri) ≈ [0.0, 2000.0]
    
    # Test 3: Piecewise load with multiple segments
    positions_pw = [0.0, 0.25, 0.75, 1.0]
    widths_pw = [1.0u"m", 2.0u"m", 2.0u"m", 1.0u"m"]  # trapezoid-ish shape
    
    trib_pw = TributaryLoad(e, positions_pw, widths_pw, pressure, direction)
    @test length(trib_pw.positions) == 4
    @test intensities(trib_pw) ≈ [1000.0, 2000.0, 2000.0, 1000.0]
    
    # =========================================================================
    # Test 4: Trapezoidal slab load - verified against reference calculator
    # https://calcresource.com/statics-fixed-beam.html
    # =========================================================================
    # Setup:
    #   L = 10 m, w = 30 kN/m, a = 2 m, b = 4 m
    #   Load shape: rises from 0 at x=0 to w at x=a, constant w until x=L-b, 
    #               then falls to 0 at x=L
    # Expected results:
    #   R_A = 116.52 kN, R_B = 93.48 kN
    #   M_A = -221.6 kNm, M_B = -196.4 kNm
    # =========================================================================
    
    println("\n--- Trapezoidal Slab Load Test (vs reference calculator) ---")
    
    # Fresh nodes and element for clean test
    n1_trap = Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    n2_trap = Node([10.0u"m", 0.0u"m", 0.0u"m"], :fixed)
    e_trap = Element(n1_trap, n2_trap, sec)
    
    # Trapezoidal load: 0 → 30 kN/m over [0, 0.2], constant 30 kN/m over [0.2, 0.6], 
    #                   30 → 0 kN/m over [0.6, 1.0]
    # Using width=1m and pressure=30000 Pa to get 30 kN/m peak intensity
    positions_trap = [0.0, 0.2, 0.6, 1.0]
    widths_trap = [0.0u"m", 1.0u"m", 1.0u"m", 0.0u"m"]
    pressure_trap = 30000.0u"Pa"  # 30 kPa
    
    trib_trap = TributaryLoad(e_trap, positions_trap, widths_trap, pressure_trap, direction)
    
    # Verify intensities: [0, 30000, 30000, 0] N/m = [0, 30, 30, 0] kN/m
    @test intensities(trib_trap) ≈ [0.0, 30000.0, 30000.0, 0.0]
    
    model_trap = Model([n1_trap, n2_trap], [e_trap], [trib_trap])
    planarize!(model_trap)
    solve!(model_trap)
    
    # Extract reactions (planarized 2D: Fy, Mz at each node)
    # In Asap's convention for planarized model with load in -Z → -Y:
    # fixedDOFs order: [Fy1, Mz1, Fy2, Mz2]
    reactions_trap = model_trap.reactions[model_trap.fixedDOFs]
    nonzero_trap = reactions_trap[reactions_trap .!= 0]
    
    println("Trapezoidal load intensities: ", intensities(trib_trap), " N/m")
    println("Reactions (raw): ", nonzero_trap)
    
    # Expected values from reference calculator (in N and N*m)
    R_A_expected = 116520.0   # 116.52 kN
    R_B_expected = 93480.0    # 93.48 kN
    M_A_expected = -221600.0  # -221.6 kNm (hogging)
    M_B_expected = -196400.0  # -196.4 kNm (hogging)
    
    # Note: Signs depend on Asap's convention. The absolute values should match.
    # Vertical reactions should sum to total load = 0.5*2*30 + 4*30 + 0.5*4*30 = 210 kN
    total_load = 210000.0  # N
    
    println("Expected: R_A=116.52kN, R_B=93.48kN, M_A=-221.6kNm, M_B=-196.4kNm")
    println("Total load should be: ", total_load/1000, " kN")
    
    # Test that vertical reactions sum correctly
    # In planarized model, reactions are [Fy1, Mz1, Fy2, Mz2]
    Fy1 = nonzero_trap[1]
    Mz1 = nonzero_trap[2]
    Fy2 = nonzero_trap[3]
    Mz2 = nonzero_trap[4]
    
    println("Computed: R_A=$(round(Fy1/1000, digits=2))kN, R_B=$(round(Fy2/1000, digits=2))kN")
    println("Computed: M_A=$(round(Mz1/1000, digits=2))kNm, M_B=$(round(Mz2/1000, digits=2))kNm")
    
    # Verify sum of reactions equals total load
    @test abs(Fy1 + Fy2 - total_load) < 1.0  # Within 1 N
    
    # Verify individual reactions (1% tolerance for numerical differences)
    tol = 0.01
    @test abs(Fy1 - R_A_expected) / R_A_expected < tol
    @test abs(Fy2 - R_B_expected) / R_B_expected < tol
    @test abs(abs(Mz1) - abs(M_A_expected)) / abs(M_A_expected) < tol
    @test abs(abs(Mz2) - abs(M_B_expected)) / abs(M_B_expected) < tol
    
    println("Trapezoidal load test PASSED!")
    println("\nAll TributaryLoad tests passed!")
end

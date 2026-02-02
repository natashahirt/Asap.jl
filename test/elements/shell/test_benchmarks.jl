#=
Shell Element Benchmark Tests
=============================

Classical benchmark problems adapted from FinEtools and FEM literature.
These tests verify the accuracy of the shell element implementations.

References:
[1] Felippa & Bergan, "A Triangular Plate Bending Element Based on an 
    Energy-Orthogonal Free Formulation", 1986, Table 3
[2] MacNeal & Harder, "A Proposed Standard Set of Problems to Test 
    Finite Element Accuracy", Finite Elements in Analysis Design, 1985
[3] Scordelis-Lo Roof - NAFEMS benchmark

Original FinEtools tests by Petr Krysl (FinEtoolsFlexStructures.jl)
Adapted for Asap.jl - MIT License
=#

using Test
using Asap
using LinearAlgebra
using SparseArrays
using Unitful

# Helper - define at top so it's available in all testsets
_test_mean(x) = sum(x) / length(x)

@testset "Shell Element Benchmarks" begin

    @testset "Simply Supported Plate - Uniform Load (Navier Solution)" begin
        # Classic plate bending benchmark: Simply supported square plate with UDL
        # Analytical solution from Felippa & Bergan 1986 (Table 3):
        #   w_center = -4.06235e-3 × p × a⁴ / D
        # where D = E×t³/(12(1-ν²))
        
        # Material and geometry
        E = 30e6    # psi (typical concrete/aluminum)
        ν = 0.3
        L = 10.0    # Square plate side length
        tL_ratio = 0.1  # thickness/length = 10%
        t = L * tL_ratio
        
        # Flexural rigidity
        D = E * t^3 / (12 * (1 - ν^2))
        
        # Load (scaled by thickness ratio for numerical stability)
        p = 1.0 * tL_ratio
        
        # Analytical center deflection (Navier solution)
        w_analytical = -4.06235e-3 * p * L^4 / D
        
        # Create mesh - 8×8 grid of triangles
        n = 8
        
        # Generate nodes
        nodes = Asap.Node[]
        for j in 0:n
            for i in 0:n
                x = (i/n - 0.5) * L  # Center at origin
                y = (j/n - 0.5) * L
                push!(nodes, Asap.Node([x*u"m", y*u"m", 0.0u"m"], :free))
            end
        end
        
        # Assign node IDs (6 DOF per node)
        for (idx, node) in enumerate(nodes)
            node.nodeID = idx
            node.globalID = collect((idx-1)*6+1 : idx*6)
        end
        
        n_nodes = length(nodes)
        n_dof = n_nodes * 6
        
        # Create triangular elements (2 triangles per cell)
        elements = Asap.ShellTri3[]
        for j in 0:n-1
            for i in 0:n-1
                # Node indices (1-based)
                n1 = j * (n+1) + i + 1
                n2 = j * (n+1) + i + 2
                n3 = (j+1) * (n+1) + i + 1
                n4 = (j+1) * (n+1) + i + 2
                
                # Triangle 1: n1-n2-n4
                tri1 = Asap.ShellTri3(
                    (nodes[n1], nodes[n2], nodes[n4]),
                    t*u"m", E*u"Pa", ν
                )
                push!(elements, tri1)
                
                # Triangle 2: n1-n4-n3
                tri2 = Asap.ShellTri3(
                    (nodes[n1], nodes[n4], nodes[n3]),
                    t*u"m", E*u"Pa", ν
                )
                push!(elements, tri2)
            end
        end
        
        # Process elements
        for elem in elements
            Asap.process!(elem)
            elem.globalID = vcat([n.globalID for n in elem.nodes]...)
        end
        
        # Assemble global stiffness
        K = zeros(n_dof, n_dof)
        for elem in elements
            gid = elem.globalID
            K[gid, gid] .+= elem.K
        end
        
        # Apply boundary conditions (soft simple support)
        # Simply supported: w = 0 on all edges
        # Constrain in-plane DOFs minimally (to prevent rigid body motion)
        fixed_dofs = Int[]
        tol = L / n / 100
        
        for node in nodes
            x = ustrip(u"m", node.position[1])
            y = ustrip(u"m", node.position[2])
            on_x_edge = abs(abs(x) - L/2) < tol
            on_y_edge = abs(abs(y) - L/2) < tol
            
            # All edges: w = 0 (deflection), θz = 0 (drilling)
            if on_x_edge || on_y_edge
                append!(fixed_dofs, node.globalID[[3, 6]])  # w, θz
            end
        end
        
        # Pin one corner to prevent rigid body in-plane motion
        corner_node = nodes[1]  # Bottom-left corner
        append!(fixed_dofs, corner_node.globalID[[1, 2]])  # u, v at corner
        
        # Constrain rotation at opposite corner to prevent spin
        opposite_corner = nodes[end]  # Top-right corner  
        push!(fixed_dofs, opposite_corner.globalID[6])  # θz (redundant but safe)
        
        fixed_dofs = unique(fixed_dofs)
        free_dofs = setdiff(1:n_dof, fixed_dofs)
        
        # Apply uniform pressure load (negative z = downward)
        F = zeros(n_dof)
        for elem in elements
            # Equivalent nodal forces: P/3 at each node
            pressure_force = -p * elem.area / 3
            for node in elem.nodes
                F[node.globalID[3]] += pressure_force  # Z-direction
            end
        end
        
        # Solve
        K_ff = K[free_dofs, free_dofs]
        F_f = F[free_dofs]
        
        u = zeros(n_dof)
        u[free_dofs] = K_ff \ F_f
        
        # Find center node displacement
        center_idx = findfirst(n -> begin
            x = ustrip(u"m", n.position[1])
            y = ustrip(u"m", n.position[2])
            abs(x) < tol && abs(y) < tol
        end, nodes)
        
        w_center = u[nodes[center_idx].globalID[3]]
        
        # Convergence check - with soft simple support, solution converges from above
        # (slightly too flexible), reaching ~100-115% at moderate mesh densities
        accuracy_percent = w_center / w_analytical * 100
        
        @test accuracy_percent > 80   # Should be close to or exceeding analytical
        @test accuracy_percent < 130  # Not overshooting too much
        @test w_center < 0  # Deflects downward
        
        @info "Simply supported plate benchmark" n=n w_analytical w_center accuracy_percent
    end
    
    @testset "Cantilever Plate - Tip Load" begin
        # Cantilever plate loaded at tip - simple bending check
        # Expected: w = P×L³/(3EI) where I = b×t³/12
        
        L = 1.0    # Length [m]
        b = 0.2    # Width [m]
        t = 0.01   # Thickness [m]
        E = 200e9  # Young's modulus [Pa]
        ν = 0.3
        P = 100.0  # Tip load [N]
        
        # Analytical beam deflection
        I = b * t^3 / 12
        w_beam = P * L^3 / (3 * E * I)
        
        # Create mesh: 4×1 elements along length
        nx = 4
        ny = 1
        
        nodes = Asap.Node[]
        for j in 0:ny
            for i in 0:nx
                x = i/nx * L
                y = j/ny * b
                push!(nodes, Asap.Node([x*u"m", y*u"m", 0.0u"m"], :free))
            end
        end
        
        for (idx, node) in enumerate(nodes)
            node.nodeID = idx
            node.globalID = collect((idx-1)*6+1 : idx*6)
        end
        
        n_nodes = length(nodes)
        n_dof = n_nodes * 6
        
        # Create elements
        elements = Asap.ShellTri3[]
        for j in 0:ny-1
            for i in 0:nx-1
                n1 = j * (nx+1) + i + 1
                n2 = j * (nx+1) + i + 2
                n3 = (j+1) * (nx+1) + i + 1
                n4 = (j+1) * (nx+1) + i + 2
                
                tri1 = Asap.ShellTri3((nodes[n1], nodes[n2], nodes[n4]), t*u"m", E*u"Pa", ν)
                tri2 = Asap.ShellTri3((nodes[n1], nodes[n4], nodes[n3]), t*u"m", E*u"Pa", ν)
                push!(elements, tri1, tri2)
            end
        end
        
        for elem in elements
            Asap.process!(elem)
            elem.globalID = vcat([n.globalID for n in elem.nodes]...)
        end
        
        # Assemble
        K = zeros(n_dof, n_dof)
        for elem in elements
            gid = elem.globalID
            K[gid, gid] .+= elem.K
        end
        
        # Fix left edge (x=0) - all DOFs
        fixed_dofs = Int[]
        tol = L / nx / 100
        for node in nodes
            x = ustrip(u"m", node.position[1])
            if abs(x) < tol
                append!(fixed_dofs, node.globalID)
            end
        end
        free_dofs = setdiff(1:n_dof, fixed_dofs)
        
        # Apply load at right edge (x=L)
        F = zeros(n_dof)
        tip_nodes = filter(n -> abs(ustrip(u"m", n.position[1]) - L) < tol, nodes)
        load_per_node = P / length(tip_nodes)
        for node in tip_nodes
            F[node.globalID[3]] = -load_per_node  # Downward
        end
        
        # Solve
        K_ff = K[free_dofs, free_dofs]
        u = zeros(n_dof)
        u[free_dofs] = K_ff \ F[free_dofs]
        
        # Get tip deflection (average of tip nodes)
        w_tip = _test_mean([u[n.globalID[3]] for n in tip_nodes])
        
        # Shell should be slightly stiffer than beam (2D vs 1D)
        # But should be within reasonable range
        @test abs(w_tip) > 0.5 * w_beam  # At least 50% of beam theory
        @test abs(w_tip) < 2.0 * w_beam  # Not more than 200%
        @test w_tip < 0  # Deflects downward
        
        @info "Cantilever plate benchmark" w_beam w_tip ratio=abs(w_tip)/w_beam
    end
    
    @testset "Patch Test - Constant Membrane Strain" begin
        # Membrane patch test: linear displacement field should produce constant strain
        # This verifies the membrane (in-plane) formulation passes basic requirements
        
        # Create 2×2 grid (4 triangles)
        nodes = [
            Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free),
            Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free),
            Asap.Node([1.0u"m", 1.0u"m", 0.0u"m"], :free),
            Asap.Node([0.0u"m", 1.0u"m", 0.0u"m"], :free),
        ]
        
        for (i, n) in enumerate(nodes)
            n.nodeID = i
            n.globalID = collect((i-1)*6+1 : i*6)
        end
        
        E = 200e9
        ν = 0.3
        t = 0.1
        
        elem1 = Asap.ShellTri3((nodes[1], nodes[2], nodes[4]), t*u"m", E*u"Pa", ν)
        elem2 = Asap.ShellTri3((nodes[2], nodes[3], nodes[4]), t*u"m", E*u"Pa", ν)
        
        Asap.process!(elem1)
        Asap.process!(elem2)
        elem1.globalID = vcat([n.globalID for n in elem1.nodes]...)
        elem2.globalID = vcat([n.globalID for n in elem2.nodes]...)
        
        # Apply linear displacement: u_x = εxx × x
        εxx = 0.001
        u = zeros(24)
        x_coords = [0.0, 1.0, 1.0, 0.0]
        for (i, x) in enumerate(x_coords)
            u[(i-1)*6 + 1] = εxx * x
        end
        
        # Get stresses (note: these are in LOCAL element coordinates)
        σ1 = Asap.stress(elem1, u)
        σ2 = Asap.stress(elem2, u)
        
        # For constant strain field, both elements should produce same stress STATE
        # but in different local coordinates. Check principal stresses instead.
        # For uniaxial strain εxx, we expect σxx = E/(1-ν²)×εxx, σyy = ν×σxx
        
        # Sort to get principal-like ordering
        σ1_sorted = sort(σ1[1:2], rev=true)
        σ2_sorted = sort(σ2[1:2], rev=true)
        
        # Both triangles should have same principal stresses
        @test isapprox(σ1_sorted[1], σ2_sorted[1], rtol=0.1)
        @test isapprox(σ1_sorted[2], σ2_sorted[2], rtol=0.1)
        @test isapprox(σ1[3], 0.0, atol=1e6)  # τxy ≈ 0
        @test isapprox(σ2[3], 0.0, atol=1e6)  # τxy ≈ 0
        
        @info "Membrane patch test" σ1 σ2
    end
    
    @testset "Pure Bending - Constant Curvature" begin
        # Apply constant curvature, check that bending moments are recovered correctly
        
        nodes = [
            Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :free),
            Asap.Node([1.0u"m", 0.0u"m", 0.0u"m"], :free),
            Asap.Node([0.5u"m", 0.866u"m", 0.0u"m"], :free),  # Equilateral
        ]
        
        for (i, n) in enumerate(nodes)
            n.nodeID = i
            n.globalID = collect((i-1)*6+1 : i*6)
        end
        
        E = 200e9
        ν = 0.3
        t = 0.1
        
        elem = Asap.ShellTri3((nodes[1], nodes[2], nodes[3]), t*u"m", E*u"Pa", ν)
        Asap.process!(elem)
        elem.globalID = vcat([n.globalID for n in elem.nodes]...)
        
        # Apply constant curvature about x-axis: κxx = 0.001 1/m
        # θy = κxx × x (rotation about y increases with x)
        κxx = 0.001
        u = zeros(18)
        x_coords = [0.0, 1.0, 0.5]
        for (i, x) in enumerate(x_coords)
            u[(i-1)*6 + 5] = κxx * x  # θy DOF
        end
        
        M = Asap.bending_moments(elem, u)
        
        # Expected: Mxx = D × κxx
        D = E * t^3 / (12 * (1 - ν^2))
        Mxx_expected = D * κxx
        
        # Check Mxx is close to expected (allow 20% tolerance for element formulation)
        @test isapprox(M[1], Mxx_expected, rtol=0.2)
        
        @info "Pure bending test" M Mxx_expected
    end
    
    @testset "Cook's Membrane - Classic Benchmark" begin
        # Cook's trapezoidal membrane - a classic test for membrane elements
        # Reference: Cook et al., "Concepts and Applications of FEA", 4th Ed.
        # Analytical tip deflection: 23.97 (converged value)
        
        E = 1.0      # Young's modulus
        ν = 1.0/3    # Poisson's ratio
        t = 1.0      # Thickness
        
        # Geometry: Trapezoidal membrane
        # Left edge (x=0): height 44, clamped
        # Right edge (x=48): height 16, loaded with shear
        width = 48.0
        height = 44.0
        free_height = 16.0
        
        # Reference tip deflection (converged value)
        u_tip_ref = 23.97
        
        # Applied load
        magn = 1.0 / free_height  # Total force = 1.0
        
        # Create mesh (16×16 - balances accuracy (~90%) with test time)
        n = 16
        
        # Generate distorted trapezoidal mesh (matching FinEtools geometry exactly)
        # The transformation creates a tilted trapezoid where:
        # - Left edge: x=0, y ranges from 0 to 44
        # - Right edge: x=48, y ranges from 44 to 60 (16 units high)
        nodes = Asap.Node[]
        for j in 0:n
            for i in 0:n
                x = i / n * width
                y_orig = j / n * height  # Original rectangular grid y
                # FinEtools transformation: y = y_orig + (x/width)*(height - y_orig/height*(height-free_height))
                y = y_orig + (x/width) * (height - y_orig/height * (height - free_height))
                push!(nodes, Asap.Node([x*u"m", y*u"m", 0.0u"m"], :free))
            end
        end
        
        for (idx, node) in enumerate(nodes)
            node.nodeID = idx
            node.globalID = collect((idx-1)*6+1 : idx*6)
        end
        
        n_nodes = length(nodes)
        n_dof = n_nodes * 6
        
        # Create triangular elements
        elements = Asap.ShellTri3[]
        for j in 0:n-1
            for i in 0:n-1
                n1 = j * (n+1) + i + 1
                n2 = j * (n+1) + i + 2
                n3 = (j+1) * (n+1) + i + 1
                n4 = (j+1) * (n+1) + i + 2
                
                tri1 = Asap.ShellTri3((nodes[n1], nodes[n2], nodes[n4]), t*u"m", E*u"Pa", ν)
                tri2 = Asap.ShellTri3((nodes[n1], nodes[n4], nodes[n3]), t*u"m", E*u"Pa", ν)
                push!(elements, tri1, tri2)
            end
        end
        
        for elem in elements
            Asap.process!(elem)
            elem.globalID = vcat([n.globalID for n in elem.nodes]...)
        end
        
        # Assemble
        K = zeros(n_dof, n_dof)
        for elem in elements
            gid = elem.globalID
            K[gid, gid] .+= elem.K
        end
        
        # Fix left edge (x=0) - clamped
        fixed_dofs = Int[]
        tol = width / n / 100
        for node in nodes
            x = ustrip(u"m", node.position[1])
            if abs(x) < tol
                append!(fixed_dofs, node.globalID)  # All 6 DOFs
            end
        end
        free_dofs = setdiff(1:n_dof, fixed_dofs)
        
        # Apply shear load on right edge (x = width)
        # Use consistent nodal forces based on tributary length
        F = zeros(n_dof)
        
        # Get right edge nodes sorted by y-coordinate
        right_edge_nodes = filter(nd -> abs(ustrip(u"m", nd.position[1]) - width) < tol, nodes)
        sort!(right_edge_nodes, by = nd -> ustrip(u"m", nd.position[2]))
        
        # Compute tributary length for each node (trapezoidal rule)
        for (i, node) in enumerate(right_edge_nodes)
            if i == 1
                # Bottom node: half distance to next node
                y_next = ustrip(u"m", right_edge_nodes[2].position[2])
                y_curr = ustrip(u"m", node.position[2])
                trib_length = (y_next - y_curr) / 2
            elseif i == length(right_edge_nodes)
                # Top node: half distance from previous node
                y_prev = ustrip(u"m", right_edge_nodes[end-1].position[2])
                y_curr = ustrip(u"m", node.position[2])
                trib_length = (y_curr - y_prev) / 2
            else
                # Interior node: half distance to neighbors
                y_prev = ustrip(u"m", right_edge_nodes[i-1].position[2])
                y_next = ustrip(u"m", right_edge_nodes[i+1].position[2])
                trib_length = (y_next - y_prev) / 2
            end
            F[node.globalID[2]] = magn * trib_length  # Y-direction shear
        end
        
        # Solve
        K_ff = K[free_dofs, free_dofs]
        u = zeros(n_dof)
        u[free_dofs] = K_ff \ F[free_dofs]
        
        # Find tip deflection at mid-height of right edge
        # Right edge goes from y=44 to y=60, so mid-point is at y=52
        mid_y = 52.0
        tip_node = nothing
        min_dist = Inf
        for node in nodes
            x = ustrip(u"m", node.position[1])
            y = ustrip(u"m", node.position[2])
            if abs(x - width) < tol
                dist = abs(y - mid_y)
                if dist < min_dist
                    min_dist = dist
                    tip_node = node
                end
            end
        end
        
        u_tip = u[tip_node.globalID[2]]  # Y-displacement
        accuracy_percent = u_tip / u_tip_ref * 100
        
        # At n=16 with correct geometry, expect ~90% accuracy (FinEtools gets ~95%)
        @test accuracy_percent > 80   # At least 80%
        @test accuracy_percent < 110  # Not overshooting
        
        @info "Cook's membrane benchmark" n=n u_tip_ref u_tip accuracy_percent
    end
    
    @testset "Composite Plate - Symmetric Layup [0/90/90/0]" begin
        # Composite plate with symmetric [0/90/90/0] layup under uniform load
        # Tests the CompositeShellTri3 element with laminate material
        
        # E-glass/epoxy properties (from FinEtools example)
        E1 = 60000.0  # Fiber direction [psi]
        E2 = 15000.0  # Transverse [psi]
        G12 = 6200.0  # In-plane shear [psi]
        G13 = G12
        G23 = 5000.0
        ν12 = 0.28
        
        # Geometry
        L = 12.0       # Side length [in]
        total_t = 0.2  # Total thickness [in]
        ply_t = total_t / 4  # 4 plies
        
        # Create laminate [0/90/90/0] - symmetric
        ply_0 = Asap.Ply("E-glass_0", E1, E2, G12, ν12, ply_t, 0.0; G13=G13, G23=G23)
        ply_90 = Asap.Ply("E-glass_90", E1, E2, G12, ν12, ply_t, 90.0; G13=G13, G23=G23)
        laminate = Asap.Laminate("0/90/90/0", [ply_0, ply_90, ply_90, ply_0])
        
        # Verify symmetric laminate has B ≈ 0
        A, B, D = Asap.laminate_stiffnesses(laminate)
        @test maximum(abs.(B)) < 1e-6 * maximum(abs.(A))  # B matrix should be ~0
        
        # Create mesh (6×6 for faster test)
        n = 6
        nodes = Asap.Node[]
        for j in 0:n
            for i in 0:n
                x = (i/n - 0.5) * L
                y = (j/n - 0.5) * L
                push!(nodes, Asap.Node([x*u"m", y*u"m", 0.0u"m"], :free))
            end
        end
        
        for (idx, node) in enumerate(nodes)
            node.nodeID = idx
            node.globalID = collect((idx-1)*6+1 : idx*6)
        end
        
        n_nodes = length(nodes)
        n_dof = n_nodes * 6
        
        # Create composite shell elements
        elements = Asap.CompositeShellTri3[]
        for j in 0:n-1
            for i in 0:n-1
                n1 = j * (n+1) + i + 1
                n2 = j * (n+1) + i + 2
                n3 = (j+1) * (n+1) + i + 1
                n4 = (j+1) * (n+1) + i + 2
                
                tri1 = Asap.CompositeShellTri3((nodes[n1], nodes[n2], nodes[n4]), laminate)
                tri2 = Asap.CompositeShellTri3((nodes[n1], nodes[n4], nodes[n3]), laminate)
                push!(elements, tri1, tri2)
            end
        end
        
        for elem in elements
            Asap.process!(elem)
            elem.globalID = vcat([n.globalID for n in elem.nodes]...)
        end
        
        # Assemble
        K = zeros(n_dof, n_dof)
        for elem in elements
            gid = elem.globalID
            K[gid, gid] .+= elem.K
        end
        
        # Simply supported edges (soft support)
        fixed_dofs = Int[]
        tol = L / n / 100
        for node in nodes
            x = ustrip(u"m", node.position[1])
            y = ustrip(u"m", node.position[2])
            on_edge = abs(abs(x) - L/2) < tol || abs(abs(y) - L/2) < tol
            if on_edge
                append!(fixed_dofs, node.globalID[[3, 6]])  # w, θz
            end
        end
        # Pin corner
        append!(fixed_dofs, nodes[1].globalID[[1, 2]])
        fixed_dofs = unique(fixed_dofs)
        free_dofs = setdiff(1:n_dof, fixed_dofs)
        
        # Uniform pressure load
        q = 0.05  # pressure [psi]
        F = zeros(n_dof)
        for elem in elements
            pressure_force = -q * elem.area / 3
            for node in elem.nodes
                F[node.globalID[3]] += pressure_force
            end
        end
        
        # Solve
        K_ff = K[free_dofs, free_dofs]
        u = zeros(n_dof)
        u[free_dofs] = K_ff \ F[free_dofs]
        
        # Get center deflection
        center_idx = findfirst(nd -> begin
            x = ustrip(u"m", nd.position[1])
            y = ustrip(u"m", nd.position[2])
            abs(x) < tol && abs(y) < tol
        end, nodes)
        
        w_center = u[nodes[center_idx].globalID[3]]
        
        # Should deflect downward and have reasonable magnitude
        @test w_center < 0
        @test abs(w_center) > 0  # Non-zero deflection
        
        @info "Composite plate benchmark" w_center laminate_thickness=Asap.thickness(laminate)
    end

end

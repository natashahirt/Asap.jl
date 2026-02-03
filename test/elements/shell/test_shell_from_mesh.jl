#=
Test Shell() creation from meshes
=================================

Tests the new Shell(mesh, section, supports) API for both:
1. DelaunayTriangulation.Triangulation (2D meshes)
2. Meshes.SimpleMesh (3D meshes)

Includes barrel vault benchmark to verify structural analysis works correctly.
=#

using Test
using Asap
using Unitful
using LinearAlgebra
import DelaunayTriangulation as DT
import Meshes

# ============================================================================
# Helper: Create barrel vault mesh data
# ============================================================================

"""
Create barrel vault mesh as (points, triangles) arrays.

Geometry (Scordelis-Lo benchmark):
- R: Radius
- L: Length (uses L/2 with symmetry)
- θ_max: Arc angle in radians
- n: Mesh density
"""
function barrel_vault_mesh_data(; R=25.0, L=50.0, θ_max=40/360*2π, n=4)
    # Create grid of points on cylindrical surface
    points = NTuple{3, Float64}[]
    node_grid = Matrix{Int}(undef, n+1, n+1)
    
    idx = 0
    for j in 0:n
        for i in 0:n
            idx += 1
            θ = i/n * θ_max
            y = j/n * (L/2)
            x = R * sin(θ)
            z = R * (cos(θ) - 1)  # z=0 at apex
            push!(points, (x, y, z))
            node_grid[i+1, j+1] = idx
        end
    end
    
    # Create triangles (2 per grid cell)
    triangles = NTuple{3, Int}[]
    for j in 1:n
        for i in 1:n
            n1 = node_grid[i, j]
            n2 = node_grid[i+1, j]
            n3 = node_grid[i, j+1]
            n4 = node_grid[i+1, j+1]
            push!(triangles, (n1, n2, n4))
            push!(triangles, (n1, n4, n3))
        end
    end
    
    # Identify boundary nodes for BCs
    diaphragm_nodes = Int[]  # y ≈ 0
    symmetry_y_nodes = Int[]  # y ≈ L/2
    apex_nodes = Int[]  # x ≈ 0
    
    tol = R / n / 100
    for (i, pt) in enumerate(points)
        x, y, z = pt
        if abs(y) < tol
            push!(diaphragm_nodes, i)
        end
        if abs(y - L/2) < tol
            push!(symmetry_y_nodes, i)
        end
        if abs(x) < tol
            push!(apex_nodes, i)
        end
    end
    
    return (
        points = points,
        triangles = triangles,
        n = n,
        R = R,
        L = L,
        θ_max = θ_max,
        diaphragm_nodes = diaphragm_nodes,
        symmetry_y_nodes = symmetry_y_nodes,
        apex_nodes = apex_nodes,
        node_grid = node_grid
    )
end

# ============================================================================
# Test: Shell from raw points + triangles
# ============================================================================

@testset "Shell from raw points + triangles" begin
    
    @testset "Basic creation from tuples" begin
        # Simple flat mesh
        points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1.0, 0.0), (1.5, 1.0, 0.0)]
        triangles = [(1, 2, 3), (2, 4, 3)]
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        # Default: :xy_plane pins all nodes at z≈0
        shells = Shell(points, triangles, section)
        @test length(shells) == 2
        
        nodes = get_nodes(shells)
        @test length(nodes) == 4
        
        # All nodes should be pinned (z=0)
        @test all(n -> n.dof == [false, false, false, true, true, true], nodes)
    end
    
    @testset "Explicit support indices" begin
        points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1.0, 0.0)]
        triangles = [(1, 2, 3)]
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        # Only pin node 1
        shells = Shell(points, triangles, section, [1])
        nodes = get_nodes(shells)
        
        @test nodes[1].dof == [false, false, false, true, true, true]  # pinned
        @test nodes[2].dof == [true, true, true, true, true, true]     # free
        @test nodes[3].dof == [true, true, true, true, true, true]     # free
    end
    
    @testset "Support type :none" begin
        points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1.0, 0.0)]
        triangles = [(1, 2, 3)]
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        shells = Shell(points, triangles, section, :none)
        nodes = get_nodes(shells)
        
        # All nodes should be free
        @test all(n -> all(n.dof), nodes)
    end
    
    @testset "From matrix" begin
        # Points as 4×3 matrix
        points = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.5 1.0 0.0; 1.5 1.0 0.0]
        triangles = [(1, 2, 3), (2, 4, 3)]
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        shells = Shell(points, triangles, section, [1, 2])
        @test length(shells) == 2
        
        nodes = get_nodes(shells)
        @test nodes[1].dof == [false, false, false, true, true, true]  # pinned
        @test nodes[2].dof == [false, false, false, true, true, true]  # pinned
        @test nodes[3].dof == [true, true, true, true, true, true]     # free
    end
    
    @testset "From vectors of vectors" begin
        points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]]
        triangles = [[1, 2, 3]]
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        shells = Shell(points, triangles, section)
        @test length(shells) == 1
        @test length(get_nodes(shells)) == 3
    end
    
    @testset "With unitful coordinates" begin
        points = [(0.0u"m", 0.0u"m", 0.0u"m"), (1.0u"m", 0.0u"m", 0.0u"m"), (0.5u"m", 1.0u"m", 0.0u"m")]
        triangles = [(1, 2, 3)]
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        shells = Shell(points, triangles, section)
        @test length(shells) == 1
        
        # Verify coordinates were converted correctly
        nodes = get_nodes(shells)
        @test ustrip(u"m", nodes[1].position[1]) ≈ 0.0
        @test ustrip(u"m", nodes[2].position[1]) ≈ 1.0
    end
    
    @testset "Plane support types" begin
        # Mesh with nodes at different z values
        points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 0.0, 1.0)]
        triangles = [(1, 2, 3)]
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        # :xy_plane should only pin z≈0
        shells = Shell(points, triangles, section, :xy_plane)
        nodes = get_nodes(shells)
        
        @test nodes[1].dof == [false, false, false, true, true, true]  # pinned (z=0)
        @test nodes[2].dof == [false, false, false, true, true, true]  # pinned (z=0)
        @test nodes[3].dof == [true, true, true, true, true, true]     # free (z=1)
        
        # :xz_plane should pin y≈0 (all nodes in this case)
        shells2 = Shell(points, triangles, section, :xz_plane)
        nodes2 = get_nodes(shells2)
        @test all(n -> n.dof == [false, false, false, true, true, true], nodes2)
    end
end

# ============================================================================
# Test: Shell from Meshes.SimpleMesh
# ============================================================================

@testset "Shell from Meshes.SimpleMesh" begin
    
    @testset "Basic creation" begin
        # Simple flat triangular mesh
        pts = [Meshes.Point(0.0, 0.0, 0.0),
               Meshes.Point(1.0, 0.0, 0.0),
               Meshes.Point(0.5, 1.0, 0.0),
               Meshes.Point(1.5, 1.0, 0.0)]
        
        conn = [Meshes.connect((1, 2, 3)), Meshes.connect((2, 4, 3))]
        mesh = Meshes.SimpleMesh(pts, conn)
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        # Default: :xy_plane pins all nodes at z≈0
        shells = Shell(mesh, section)
        @test length(shells) == 2
        
        nodes = get_nodes(shells)
        @test length(nodes) == 4
        
        # All nodes should be pinned (z=0)
        @test all(n -> n.dof == [false, false, false, true, true, true], nodes)
    end
    
    @testset "Explicit support indices" begin
        pts = [Meshes.Point(0.0, 0.0, 0.0),
               Meshes.Point(1.0, 0.0, 0.0),
               Meshes.Point(0.5, 1.0, 0.0)]
        
        conn = [Meshes.connect((1, 2, 3))]
        mesh = Meshes.SimpleMesh(pts, conn)
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        # Only pin node 1
        shells = Shell(mesh, section, [1])
        nodes = get_nodes(shells)
        
        @test nodes[1].dof == [false, false, false, true, true, true]  # pinned
        @test nodes[2].dof == [true, true, true, true, true, true]     # free
        @test nodes[3].dof == [true, true, true, true, true, true]     # free
    end
    
    @testset "Support type :none" begin
        pts = [Meshes.Point(0.0, 0.0, 0.0),
               Meshes.Point(1.0, 0.0, 0.0),
               Meshes.Point(0.5, 1.0, 0.0)]
        
        conn = [Meshes.connect((1, 2, 3))]
        mesh = Meshes.SimpleMesh(pts, conn)
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        shells = Shell(mesh, section, :none)
        nodes = get_nodes(shells)
        
        # All nodes should be free
        @test all(n -> all(n.dof), nodes)
    end
    
    @testset "Plane support types" begin
        # Mesh with nodes at different z values
        pts = [Meshes.Point(0.0, 0.0, 0.0),   # z=0
               Meshes.Point(1.0, 0.0, 0.0),   # z=0
               Meshes.Point(0.5, 0.0, 1.0)]   # z=1
        
        conn = [Meshes.connect((1, 2, 3))]
        mesh = Meshes.SimpleMesh(pts, conn)
        
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        # :xy_plane should only pin z≈0
        shells = Shell(mesh, section, :xy_plane)
        nodes = get_nodes(shells)
        
        @test nodes[1].dof == [false, false, false, true, true, true]  # pinned (z=0)
        @test nodes[2].dof == [false, false, false, true, true, true]  # pinned (z=0)
        @test nodes[3].dof == [true, true, true, true, true, true]     # free (z=1)
    end
end

# ============================================================================
# Test: Shell from DT.Triangulation
# ============================================================================

@testset "Shell from DT.Triangulation" begin
    
    @testset "Basic creation with supports" begin
        # Create corner nodes for triangulation
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([1.0u"m", 1.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 1.0u"m", 0.0u"m"], :pinned)
        
        tri = mesh((n1, n2, n3, n4), 2)
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        # Pin corners only
        shells = Shell(tri, section, [1, 2, 3, 4])
        nodes = get_nodes(shells)
        
        @test length(shells) > 0
        @test length(nodes) > 4  # corners + interior nodes
        
        # First 4 nodes (corners) should be pinned
        pinned = filter(n -> n.dof == [false, false, false, true, true, true], nodes)
        @test length(pinned) == 4
    end
    
    @testset ":xy_plane pins all at z=0" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Node([0.5u"m", 1.0u"m", 0.0u"m"], :free)
        
        tri = mesh((n1, n2, n3), 2)
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        # :xy_plane at z=0 should pin all nodes
        shells = Shell(tri, section, :xy_plane; z=0.0)
        nodes = get_nodes(shells)
        
        @test all(n -> n.dof == [false, false, false, true, true, true], nodes)
    end
    
    @testset ":none creates all free" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Node([1.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Node([0.5u"m", 1.0u"m", 0.0u"m"], :free)
        
        tri = mesh((n1, n2, n3), 2)
        section = ShellSection(0.1u"m", 30e9u"Pa", 0.2)
        
        shells = Shell(tri, section, :none)
        nodes = get_nodes(shells)
        
        @test all(n -> all(n.dof), nodes)
    end
end

# ============================================================================
# Test: Barrel Vault with Meshes.SimpleMesh
# ============================================================================

@testset "Barrel Vault - Meshes.SimpleMesh workflow" begin
    data = barrel_vault_mesh_data(n=4)
    
    # Create Meshes.jl SimpleMesh
    pts = [Meshes.Point(p...) for p in data.points]
    conn = [Meshes.connect(t) for t in data.triangles]
    mesh = Meshes.SimpleMesh(pts, conn)
    
    # Section (Scordelis-Lo parameters)
    # Note: Original benchmark uses ν=0.0, but ShellSection requires ν > 0
    # Using ν=0.001 which is effectively zero
    thickness = 0.25  # ft
    E = 4.32e8        # psf
    ν = 0.001
    section = ShellSection(thickness*u"m", E*u"Pa", ν)
    
    # Create shells with diaphragm supports (y≈0)
    shells = Shell(mesh, section, data.diaphragm_nodes)
    nodes = get_nodes(shells)
    
    @test length(shells) == 2 * data.n^2  # 2 triangles per cell
    @test length(nodes) == (data.n + 1)^2
    
    # Note: nodes from get_nodes() are in order they appear in shells (mesh vertex order)
    # so nodes[i] corresponds to mesh vertex i
    
    # Verify diaphragm nodes are pinned (dof = [false, false, false, true, true, true])
    pinned_count = count(n -> n.dof == [false, false, false, true, true, true], nodes)
    @test pinned_count == length(data.diaphragm_nodes)
    
    # Apply additional BCs manually for full Scordelis-Lo
    # Use y-position to identify symmetry nodes
    tol = data.R / data.n / 100
    for node in nodes
        y = ustrip(u"m", node.position[2])
        x = ustrip(u"m", node.position[1])
        
        # Symmetry at y=L/2: fix y, θx, θz
        if abs(y - data.L/2) < tol
            node.dof[[2, 4, 6]] .= false
        end
        # Apex symmetry at x=0: fix x, θy, θz
        if abs(x) < tol
            node.dof[[1, 5, 6]] .= false
        end
    end
    
    # Create load (uniform pressure) - AreaLoad takes scalar pressure + direction
    p = 90.0  # psf magnitude
    loads = [AreaLoad(shells, p*u"Pa"; direction=(0.0, 0.0, -1.0))]
    
    # Build and solve model
    model = ShellModel(nodes, shells, loads)
    solve!(model)
    
    # Check that we got a solution
    @test model.u !== nothing
    @test !all(iszero, model.u)
    
    # Find deflection at free edge by position (i=n, j=n corresponds to x=max, y=L/2)
    # The free edge is at θ=θ_max (max x), y=L/2
    x_max = data.R * sin(data.θ_max)
    free_edge_node = nothing
    for node in nodes
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        if abs(x - x_max) < tol && abs(y - data.L/2) < tol
            free_edge_node = node
            break
        end
    end
    
    @test free_edge_node !== nothing
    w_z = ustrip(u"m", free_edge_node.displacement[3])
    
    @info "Barrel vault (Meshes.jl): w_z = $w_z m at free edge"
    
    # Analytical solution: -0.3024 ft (negative = downward)
    # We expect ~65-70% of analytical for n=4
    @test w_z < 0  # Should deflect downward
end

# ============================================================================
# Test: Barrel Vault with manual node creation (traditional workflow)
# ============================================================================

@testset "Barrel Vault - Manual node workflow (comparison)" begin
    data = barrel_vault_mesh_data(n=4)
    
    # Traditional approach: create nodes manually with fixity
    nodes = Node[]
    tol = data.R / data.n / 100
    
    for (i, pt) in enumerate(data.points)
        x, y, z = pt
        pos = [x, y, z] .* u"m"
        
        # Determine fixity based on location
        dof = [true, true, true, true, true, true]  # start free
        
        # Rigid diaphragm at y=0
        if abs(y) < tol
            dof[[1, 3, 5]] .= false  # fix x, z, θy
        end
        # Symmetry at y=L/2
        if abs(y - data.L/2) < tol
            dof[[2, 4, 6]] .= false  # fix y, θx, θz
        end
        # Apex symmetry at x=0
        if abs(x) < tol
            dof[[1, 5, 6]] .= false  # fix x, θy, θz
        end
        
        node = Node(pos, :free)
        node.dof .= dof
        push!(nodes, node)
    end
    
    # Create shells directly
    section = ShellSection(0.25u"m", 4.32e8u"Pa", 0.001)
    shells = ShellTri3[]
    
    for (i, j, k) in data.triangles
        push!(shells, ShellTri3(
            (nodes[i], nodes[j], nodes[k]),
            section.thickness * u"m",
            section.E * u"Pa",
            section.ν;
            ρ = section.ρ
        ))
    end
    
    # Load and solve (AreaLoad takes scalar pressure + direction)
    loads = [AreaLoad(shells, 90.0u"Pa"; direction=(0.0, 0.0, -1.0))]
    model = ShellModel(nodes, shells, loads)
    solve!(model)
    
    # Get result - find free edge node by position
    x_max = data.R * sin(data.θ_max)
    free_edge_node = nothing
    for node in nodes
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        if abs(x - x_max) < tol && abs(y - data.L/2) < tol
            free_edge_node = node
            break
        end
    end
    @test free_edge_node !== nothing
    
    w_z_manual = ustrip(u"m", free_edge_node.displacement[3])
    
    @info "Barrel vault (manual): w_z = $w_z_manual m at free edge"
    
    @test w_z_manual < 0
end

# ============================================================================
# Test: Compare Meshes.SimpleMesh vs manual (should be identical)
# ============================================================================

# ============================================================================
# Test: Barrel Vault with raw points + triangles
# ============================================================================

@testset "Barrel Vault - Raw points/triangles workflow" begin
    data = barrel_vault_mesh_data(n=4)
    
    # Use raw points and triangles directly from mesh data
    section = ShellSection(0.25u"m", 4.32e8u"Pa", 0.001)
    
    # Create shells with :none support initially, then apply BCs
    shells = Shell(data.points, data.triangles, section, :none)
    nodes = get_nodes(shells)
    
    @test length(shells) == 2 * data.n^2
    @test length(nodes) == (data.n + 1)^2
    
    # Apply full Scordelis-Lo BCs
    tol = data.R / data.n / 100
    for node in nodes
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        
        # Rigid diaphragm at y=0
        if abs(y) < tol
            node.dof[[1, 3, 5]] .= false
        end
        # Symmetry at y=L/2
        if abs(y - data.L/2) < tol
            node.dof[[2, 4, 6]] .= false
        end
        # Apex symmetry at x=0
        if abs(x) < tol
            node.dof[[1, 5, 6]] .= false
        end
    end
    
    # Load and solve
    loads = [AreaLoad(shells, 90.0u"Pa"; direction=(0.0, 0.0, -1.0))]
    model = ShellModel(nodes, shells, loads)
    solve!(model)
    
    @test model.u !== nothing
    @test !all(iszero, model.u)
    
    # Find free edge deflection
    x_max = data.R * sin(data.θ_max)
    free_edge_node = nothing
    for node in nodes
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        if abs(x - x_max) < tol && abs(y - data.L/2) < tol
            free_edge_node = node
            break
        end
    end
    
    @test free_edge_node !== nothing
    w_z = ustrip(u"m", free_edge_node.displacement[3])
    @info "Barrel vault (raw points): w_z = $w_z m at free edge"
    
    @test w_z < 0  # Should deflect downward
end

# ============================================================================
# Test: Compare all three interfaces
# ============================================================================

@testset "Meshes.SimpleMesh vs manual comparison" begin
    data = barrel_vault_mesh_data(n=4)
    tol = data.R / data.n / 100
    
    # --- Meshes.jl approach ---
    pts = [Meshes.Point(p...) for p in data.points]
    conn = [Meshes.connect(t) for t in data.triangles]
    simplemesh = Meshes.SimpleMesh(pts, conn)
    
    section = ShellSection(0.25u"m", 4.32e8u"Pa", 0.001)
    # Start with :none supports, apply BCs identically to manual
    shells_mesh = Shell(simplemesh, section, :none)
    nodes_mesh = get_nodes(shells_mesh)
    
    # Apply IDENTICAL BCs as manual approach
    for node in nodes_mesh
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        
        # Rigid diaphragm at y=0
        if abs(y) < tol
            node.dof[[1, 3, 5]] .= false
        end
        # Symmetry at y=L/2
        if abs(y - data.L/2) < tol
            node.dof[[2, 4, 6]] .= false
        end
        # Apex symmetry at x=0
        if abs(x) < tol
            node.dof[[1, 5, 6]] .= false
        end
    end
    
    loads_mesh = [AreaLoad(shells_mesh, 90.0u"Pa"; direction=(0.0, 0.0, -1.0))]
    model_mesh = ShellModel(nodes_mesh, shells_mesh, loads_mesh)
    solve!(model_mesh)
    
    # --- Manual approach ---
    nodes_manual = Node[]
    for (i, pt) in enumerate(data.points)
        x, y, z = pt
        pos = [x, y, z] .* u"m"
        dof = [true, true, true, true, true, true]
        if abs(y) < tol
            dof[[1, 3, 5]] .= false
        end
        if abs(y - data.L/2) < tol
            dof[[2, 4, 6]] .= false
        end
        if abs(x) < tol
            dof[[1, 5, 6]] .= false
        end
        node = Node(pos, :free)
        node.dof .= dof
        push!(nodes_manual, node)
    end
    
    shells_manual = ShellTri3[]
    for (i, j, k) in data.triangles
        push!(shells_manual, ShellTri3(
            (nodes_manual[i], nodes_manual[j], nodes_manual[k]),
            section.thickness * u"m", section.E * u"Pa", section.ν; ρ = section.ρ
        ))
    end
    
    loads_manual = [AreaLoad(shells_manual, 90.0u"Pa"; direction=(0.0, 0.0, -1.0))]
    model_manual = ShellModel(nodes_manual, shells_manual, loads_manual)
    solve!(model_manual)
    
    # --- Compare results ---
    # Find free edge nodes by position
    x_max = data.R * sin(data.θ_max)
    
    free_edge_mesh = nothing
    for node in nodes_mesh
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        if abs(x - x_max) < tol && abs(y - data.L/2) < tol
            free_edge_mesh = node
            break
        end
    end
    
    free_edge_manual = nothing
    for node in nodes_manual
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        if abs(x - x_max) < tol && abs(y - data.L/2) < tol
            free_edge_manual = node
            break
        end
    end
    
    @test free_edge_mesh !== nothing
    @test free_edge_manual !== nothing
    
    w_mesh = ustrip(u"m", free_edge_mesh.displacement[3])
    w_manual = ustrip(u"m", free_edge_manual.displacement[3])
    
    @info "Comparison: Meshes.jl=$w_mesh, Manual=$w_manual"
    
    # Should be identical (or very close)
    @test isapprox(w_mesh, w_manual, rtol=1e-10)
end

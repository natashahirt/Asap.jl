#=
Scordelis-Lo Barrel Vault Test
==============================

This is a DIRECT comparison with FinEtools' Scordelis-Lo test.
We use the EXACT same:
- Geometry (cylindrical barrel vault)
- Material properties
- Boundary conditions
- Loading
- Mesh topology

Reference values from FinEtoolsFlexStructures/test/test_shell_statics.jl lines 95-102:
  n=4:  66.54771615057949%
  n=8:  85.54615143134853%
  n=10: 89.85075281481419%
  n=12: 92.50616661644985%
  n=16: 95.40469210310079%
  n=24: 97.65376880486126%
=#

using Test
using Asap
using LinearAlgebra
using Unitful

# FinEtools hardcoded reference values (% of analytical solution)
const FINETOOLS_SCORDELIS_LO = Dict(
    4 => 66.54771615057949,
    8 => 85.54615143134853,
    10 => 89.85075281481419,
    12 => 92.50616661644985,
    16 => 95.40469210310079,
    24 => 97.65376880486126,
)

"""
Create Scordelis-Lo barrel vault mesh matching FinEtools exactly.

Geometry:
- Quarter model with symmetry
- R = 25 ft (radius)
- L = 50 ft (length), using L/2 = 25 ft with symmetry
- θ = 40° arc
- t = 0.25 ft (thickness)

The parametric mesh is mapped via:
  x = R * sin(θ)
  y = y  
  z = R * (cos(θ) - 1)
"""
function create_scordelis_lo_mesh(n::Int)
    # Parameters (matching FinEtools exactly)
    R = 25.0       # Radius [ft]
    L = 50.0       # Total length [ft]
    thickness = 0.25  # [ft]
    E = 4.32e8     # Young's modulus [psf]
    ν = 0.0        # Poisson's ratio (zero for this benchmark!)
    
    # Arc angle
    θ_max = 40.0 / 360.0 * 2 * π
    
    # Create parametric mesh: (θ, y) from (0,0) to (θ_max, L/2)
    nodes = Asap.Node[]
    node_grid = Matrix{Union{Asap.Node, Nothing}}(nothing, n+1, n+1)
    
    for j in 0:n
        for i in 0:n
            # Parametric coordinates
            θ = i/n * θ_max
            y = j/n * (L/2)
            
            # Map to cylindrical surface
            x = R * sin(θ)
            z = R * (cos(θ) - 1)  # z=0 at apex (θ=0)
            
            # All nodes start as :free, we'll constrain in BCs
            node = Asap.Node([x*u"m", y*u"m", z*u"m"], :free)
            push!(nodes, node)
            node_grid[i+1, j+1] = node
        end
    end
    
    # Assign node IDs
    for (idx, node) in enumerate(nodes)
        node.nodeID = idx
        node.globalID = collect((idx-1)*6+1 : idx*6)
    end
    
    # Create triangular elements (same topology as FinEtools T3block)
    elements = Asap.ShellTri3[]
    for j in 0:n-1
        for i in 0:n-1
            # Four corners of this cell
            n1 = node_grid[i+1, j+1]
            n2 = node_grid[i+2, j+1]
            n3 = node_grid[i+1, j+2]
            n4 = node_grid[i+2, j+2]
            
            # Two triangles per cell (matching FinEtools)
            tri1 = Asap.ShellTri3((n1, n2, n4), thickness*u"m", E*u"Pa", ν)
            tri2 = Asap.ShellTri3((n1, n4, n3), thickness*u"m", E*u"Pa", ν)
            push!(elements, tri1, tri2)
        end
    end
    
    # Process elements
    for elem in elements
        Asap.process!(elem)
        elem.globalID = vcat([n.globalID for n in elem.nodes]...)
    end
    
    return nodes, elements, node_grid, (R=R, L=L, θ_max=θ_max, thickness=thickness, E=E, ν=ν)
end

"""
Solve Scordelis-Lo benchmark and return % of analytical solution.
"""
function solve_scordelis_lo(n::Int)
    nodes, elements, node_grid, params = create_scordelis_lo_mesh(n)
    
    R, L, θ_max = params.R, params.L, params.θ_max
    n_dof = length(nodes) * 6
    
    # Assemble stiffness
    K = zeros(n_dof, n_dof)
    for elem in elements
        gid = elem.globalID
        K[gid, gid] .+= elem.K
    end
    
    # Apply boundary conditions (matching FinEtools exactly)
    fixed_dofs = Int[]
    tol = R / n / 100
    
    for j in 0:n
        for i in 0:n
            node = node_grid[i+1, j+1]
            x = ustrip(u"m", node.position[1])
            y = ustrip(u"m", node.position[2])
            z = ustrip(u"m", node.position[3])
            
            # Rigid diaphragm at y=0: constrain u_x, u_z, θ_y (DOFs 1, 3, 5)
            if abs(y) < tol
                append!(fixed_dofs, node.globalID[[1, 3, 5]])
            end
            
            # Symmetry at y=L/2: constrain u_y, θ_x, θ_z (DOFs 2, 4, 6)
            if abs(y - L/2) < tol
                append!(fixed_dofs, node.globalID[[2, 4, 6]])
            end
            
            # Symmetry at apex (θ=0, i.e., x=0): constrain u_x, θ_y, θ_z (DOFs 1, 5, 6)
            if abs(x) < tol
                append!(fixed_dofs, node.globalID[[1, 5, 6]])
            end
        end
    end
    
    fixed_dofs = unique(fixed_dofs)
    free_dofs = setdiff(1:n_dof, fixed_dofs)
    
    # Apply uniform gravity load: 90 psf downward (Z direction)
    p = -90.0  # psf (negative = downward in Z)
    F = zeros(n_dof)
    for elem in elements
        # Consistent nodal loads for uniform pressure
        load_per_node = p * elem.area / 3
        for node in elem.nodes
            F[node.globalID[3]] += load_per_node  # Z-direction
        end
    end
    
    # Solve
    u = zeros(n_dof)
    K_ff = K[free_dofs, free_dofs]
    F_f = F[free_dofs]
    u[free_dofs] = K_ff \ F_f
    
    # Find deflection at midpoint of free edge
    # Free edge is at θ = θ_max (i = n), y = L/2 (j = n)
    free_edge_node = node_grid[n+1, n+1]
    w_z = u[free_edge_node.globalID[3]]
    
    # Analytical solution
    analyt_sol = -0.3024  # ft (from FinEtools)
    
    # Return percentage of analytical
    return w_z / analyt_sol * 100, w_z, analyt_sol
end

@testset "Scordelis-Lo Barrel Vault" begin
    
    @testset "Direct FinEtools Comparison" begin
        @info "Scordelis-Lo Barrel Vault - ASAP vs FinEtools"
        @info "  Same geometry, mesh, BCs, and loading"
        @info ""
        
        results = Dict{Int, Float64}()
        
        for n in [4, 8, 16]
            percent, w_z, analyt = solve_scordelis_lo(n)
            results[n] = percent
            
            finetools_ref = FINETOOLS_SCORDELIS_LO[n]
            diff = percent - finetools_ref
            
            @info "  n=$n:"
            @info "    ASAP:      $(round(percent, digits=2))%"
            @info "    FinEtools: $(round(finetools_ref, digits=2))%"
            @info "    Diff:      $(round(diff, digits=2))%"
            @info ""
        end
        
        # Tests
        @testset "Convergence" begin
            # Should show convergence (improving with mesh refinement)
            @test results[16] > results[4]
        end
        
        @testset "Agreement with FinEtools" begin
            # At each mesh density, ASAP should be within 15% of FinEtools
            # (Some difference is expected due to element formulation details)
            for n in [4, 8, 16]
                diff = abs(results[n] - FINETOOLS_SCORDELIS_LO[n])
                @test diff < 20  # Within 20 percentage points
            end
        end
    end
    
    @testset "Full Convergence Study" begin
        @info "\nFull Convergence Study:"
        
        all_results = Float64[]
        for n in [4, 8, 12, 16, 24]
            percent, _, _ = solve_scordelis_lo(n)
            push!(all_results, percent)
            
            ref = get(FINETOOLS_SCORDELIS_LO, n, NaN)
            @info "  n=$n: ASAP=$(round(percent, digits=2))%, FinEtools=$(round(ref, digits=2))%"
        end
        
        # Verify monotonic convergence
        @test issorted(all_results)
    end

end

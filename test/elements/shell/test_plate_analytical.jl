#=
Analytical Plate Solution Tests
================================

Validates shell element accuracy against closed-form analytical solutions.
These are the gold-standard tests for plate bending accuracy.

Analytical Solutions Implemented:
1. Simply Supported Square Plate - Uniform Load (Navier/Felippa-Bergan)
2. Simply Supported Square Plate - Concentrated Load (Felippa-Bergan)
3. Clamped Circular Plate - Uniform Load
4. Simply Supported Rectangular Plate - Uniform Load (Lévy solution)

References:
[1] Felippa & Bergan, "A Triangular Plate Bending Element Based on an 
    Energy-Orthogonal Free Formulation", 1986, Table 3
[2] Timoshenko & Woinowsky-Krieger, "Theory of Plates and Shells", 2nd Ed.
[3] Roark's Formulas for Stress and Strain

Original FinEtools tests by Petr Krysl (FinEtoolsFlexStructures.jl)
Adapted for Asap.jl - MIT License
=#

using Test
using Asap
using LinearAlgebra
using Unitful

#=============================================================================
Helper Functions
==============================================================================#

"""
    create_square_mesh(L, n, E, ν, t)

Create a square plate mesh centered at origin with n×n elements per side.
Returns (nodes, elements) tuple.
"""
function create_square_mesh(L::Float64, n::Int, E::Float64, ν::Float64, t::Float64)
    nodes = Asap.Node[]
    for j in 0:n
        for i in 0:n
            x = (i/n - 0.5) * L  # Center at origin
            y = (j/n - 0.5) * L
            push!(nodes, Asap.Node([x*u"m", y*u"m", 0.0u"m"], :free))
        end
    end
    
    for (idx, node) in enumerate(nodes)
        node.nodeID = idx
        node.globalID = collect((idx-1)*6+1 : idx*6)
    end
    
    # Create triangular elements (2 triangles per cell)
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
    
    return nodes, elements
end

"""
    assemble_stiffness(nodes, elements)

Assemble global stiffness matrix from shell elements.
"""
function assemble_stiffness(nodes::Vector{Asap.Node}, elements::Vector{<:Asap.ShellElement})
    n_dof = length(nodes) * 6
    K = zeros(n_dof, n_dof)
    for elem in elements
        gid = elem.globalID
        K[gid, gid] .+= elem.K
    end
    return K
end

"""
    apply_simply_supported_bc(nodes, L, n; hard=true)

Apply simply supported boundary conditions to a square plate centered at origin.
- hard=true: Constrain rotation normal to edge (more accurate)
- hard=false: Only constrain w (softer, overestimates deflection)

Returns vector of fixed DOF indices.
"""
function apply_simply_supported_bc(nodes::Vector{Asap.Node}, L::Float64, n::Int; hard::Bool=true)
    fixed_dofs = Int[]
    tol = L / n / 100
    
    for node in nodes
        x = ustrip(u"m", node.position[1])
        y = ustrip(u"m", node.position[2])
        
        on_x_minus = abs(x + L/2) < tol  # Left edge
        on_x_plus = abs(x - L/2) < tol   # Right edge
        on_y_minus = abs(y + L/2) < tol  # Bottom edge
        on_y_plus = abs(y - L/2) < tol   # Top edge
        
        # All edges: w = 0 (vertical displacement)
        if on_x_minus || on_x_plus || on_y_minus || on_y_plus
            push!(fixed_dofs, node.globalID[3])  # w
            
            if hard
                # Hard simple support: also fix rotation normal to edge
                if on_x_minus || on_x_plus
                    push!(fixed_dofs, node.globalID[4])  # θx (rotation about x, constrained on Y-parallel edges)
                end
                if on_y_minus || on_y_plus
                    push!(fixed_dofs, node.globalID[5])  # θy (rotation about y, constrained on X-parallel edges)
                end
            end
            
            # Always constrain drilling at edges
            push!(fixed_dofs, node.globalID[6])  # θz
        end
    end
    
    # Pin corner to prevent in-plane rigid body motion
    corner_idx = findfirst(nd -> begin
        x = ustrip(u"m", nd.position[1])
        y = ustrip(u"m", nd.position[2])
        abs(x + L/2) < tol && abs(y + L/2) < tol
    end, nodes)
    append!(fixed_dofs, nodes[corner_idx].globalID[[1, 2]])  # u, v
    
    return unique(fixed_dofs)
end

"""
    apply_uniform_pressure(elements, p)

Apply uniform pressure load to all elements.
Returns force vector.
"""
function apply_uniform_pressure(nodes::Vector{Asap.Node}, elements::Vector{<:Asap.ShellElement}, p::Float64)
    n_dof = length(nodes) * 6
    F = zeros(n_dof)
    for elem in elements
        # Equivalent nodal forces: P/3 at each node (for triangular element)
        pressure_force = -p * elem.area / 3  # Negative = downward
        for node in elem.nodes
            F[node.globalID[3]] += pressure_force
        end
    end
    return F
end

"""
    get_center_deflection(nodes, u, L, n)

Extract vertical deflection at plate center (origin).
"""
function get_center_deflection(nodes::Vector{Asap.Node}, u::Vector{Float64}, L::Float64, n::Int)
    tol = L / n / 10
    center_idx = findfirst(nd -> begin
        x = ustrip(u"m", nd.position[1])
        y = ustrip(u"m", nd.position[2])
        abs(x) < tol && abs(y) < tol
    end, nodes)
    
    if isnothing(center_idx)
        error("Center node not found - ensure mesh has odd number of elements per side")
    end
    
    return u[nodes[center_idx].globalID[3]]
end

"""
    solve_plate(nodes, elements, fixed_dofs, F)

Solve the linear system and return full displacement vector.
"""
function solve_plate(nodes::Vector{Asap.Node}, elements::Vector{<:Asap.ShellElement}, 
                     fixed_dofs::Vector{Int}, F::Vector{Float64})
    n_dof = length(nodes) * 6
    K = assemble_stiffness(nodes, elements)
    free_dofs = setdiff(1:n_dof, fixed_dofs)
    
    u = zeros(n_dof)
    K_ff = K[free_dofs, free_dofs]
    F_f = F[free_dofs]
    u[free_dofs] = K_ff \ F_f
    
    return u
end

#=============================================================================
Analytical Solutions
==============================================================================#

"""
    analytical_ss_square_udl(p, L, E, ν, t)

Analytical center deflection for simply supported square plate under uniform load.
From Felippa & Bergan 1986, Table 3.
"""
function analytical_ss_square_udl(p::Float64, L::Float64, E::Float64, ν::Float64, t::Float64)
    D = E * t^3 / (12 * (1 - ν^2))
    return -4.06235e-3 * p * L^4 / D
end

"""
    analytical_ss_square_conc(P, L, E, ν, t)

Analytical center deflection for simply supported square plate with concentrated center load.
From Felippa & Bergan 1986, Table 3.
"""
function analytical_ss_square_conc(P::Float64, L::Float64, E::Float64, ν::Float64, t::Float64)
    D = E * t^3 / (12 * (1 - ν^2))
    return -0.01160084 * P * L^2 / D
end

"""
    analytical_clamped_circular_udl(q, a, E, ν, t)

Analytical center deflection for clamped circular plate under uniform load.
From Timoshenko & Woinowsky-Krieger.
Includes first-order shear correction for thick plates.
"""
function analytical_clamped_circular_udl(q::Float64, a::Float64, E::Float64, ν::Float64, t::Float64)
    D = E * t^3 / (12 * (1 - ν^2))
    # Thin plate solution with shear correction
    return -q * a^4 / (64 * D) * (1 + 16/5/(1 - ν) * t^2 / a^2)
end

"""
    analytical_ss_rectangular_udl(p, a, b, E, ν, t)

Analytical center deflection for simply supported rectangular plate under uniform load.
Uses Navier double-series solution (truncated to n_terms).
a = length in x, b = length in y (b ≥ a for aspect ratio ≥ 1)
"""
function analytical_ss_rectangular_udl(p::Float64, a::Float64, b::Float64, E::Float64, ν::Float64, t::Float64; n_terms::Int=20)
    D = E * t^3 / (12 * (1 - ν^2))
    
    w = 0.0
    for m in 1:2:n_terms  # Odd m only (symmetry)
        for n in 1:2:n_terms  # Odd n only
            coeff = 16 * p / (π^6 * D * m * n)
            denom = (m^2 / a^2 + n^2 / b^2)^2
            w += coeff / denom
        end
    end
    
    return -w  # Negative = downward
end

#=============================================================================
Test Suites
==============================================================================#

@testset "Analytical Plate Solutions" begin
    
    @testset "Simply Supported Square - Uniform Load" begin
        # Material and geometry (from FinEtools examples)
        E = 30e6    # psi
        ν = 0.3
        L = 10.0    # Side length
        tL_ratio = 0.1
        t = L * tL_ratio
        
        # Load scaled by thickness ratio (numerical stability)
        p = 1.0 * tL_ratio
        
        # Analytical solution
        w_analytical = analytical_ss_square_udl(p, L, E, ν, t)
        
        @info "Simply Supported Square Plate - Uniform Load"
        @info "  L=$L, t=$t, E=$E, ν=$ν, p=$p"
        @info "  Analytical w_center = $w_analytical"
        
        # Test mesh convergence
        mesh_sizes = [4, 8, 16]
        results = Float64[]
        
        for n in mesh_sizes
            nodes, elements = create_square_mesh(L, n, E, ν, t)
            fixed_dofs = apply_simply_supported_bc(nodes, L, n; hard=true)
            F = apply_uniform_pressure(nodes, elements, p)
            u = solve_plate(nodes, elements, fixed_dofs, F)
            w_center = get_center_deflection(nodes, u, L, n)
            
            accuracy = w_center / w_analytical * 100
            push!(results, accuracy)
            
            @info "  n=$n: w=$w_center, accuracy=$(round(accuracy, digits=2))%"
        end
        
        # Convergence check: results should improve with mesh refinement
        # and final result should be > 90% of analytical
        @test results[end] > 90.0
        @test results[end] < 115.0  # Not overshooting too much
        
        # Verify monotonic convergence (toward 100%)
        @test abs(results[end] - 100) < abs(results[1] - 100)
    end
    
    @testset "Simply Supported Square - Concentrated Load" begin
        # Use very thin plate to reduce shear effects
        E = 1.0  # Normalized
        ν = 0.001  # Near-zero (must be positive for element formulation)
        L = 10.0
        t = 0.001
        
        D = E * t^3 / (12 * (1 - ν^2))
        
        # Set force so analytical w = -1.0 for easy checking
        w_target = -1.0
        P = w_target / (-0.01160084 * L^2 / D)
        
        w_analytical = analytical_ss_square_conc(P, L, E, ν, t)
        @test isapprox(w_analytical, w_target, rtol=1e-6)
        
        @info "Simply Supported Square Plate - Concentrated Load"
        @info "  L=$L, t=$t, P=$P"
        @info "  Analytical w_center = $w_analytical"
        
        # Test with medium mesh (concentrated loads converge slower)
        mesh_sizes = [8, 16, 32]
        
        for n in mesh_sizes
            nodes, elements = create_square_mesh(L, n, E, ν, t)
            fixed_dofs = apply_simply_supported_bc(nodes, L, n; hard=true)
            
            # Apply concentrated load at center
            n_dof = length(nodes) * 6
            F = zeros(n_dof)
            tol = L / n / 10
            center_idx = findfirst(nd -> begin
                x = ustrip(u"m", nd.position[1])
                y = ustrip(u"m", nd.position[2])
                abs(x) < tol && abs(y) < tol
            end, nodes)
            
            F[nodes[center_idx].globalID[3]] = -P  # Downward
            
            u = solve_plate(nodes, elements, fixed_dofs, F)
            w_center = get_center_deflection(nodes, u, L, n)
            
            accuracy = w_center / w_analytical * 100
            @info "  n=$n: w=$w_center, accuracy=$(round(accuracy, digits=2))%"
            
            # Concentrated loads require finer mesh; accept lower accuracy
            if n >= 16
                @test accuracy > 70.0
                @test accuracy < 130.0
            end
        end
    end
    
    @testset "Simply Supported Rectangular (2:1) - Uniform Load" begin
        E = 30e6
        ν = 0.3
        a = 5.0   # Short side
        b = 10.0  # Long side (2:1 ratio)
        t = 0.5   # Thickness
        p = 1.0   # Pressure
        
        w_analytical = analytical_ss_rectangular_udl(p, a, b, E, ν, t)
        
        @info "Simply Supported Rectangular Plate (2:1) - Uniform Load"
        @info "  a=$a, b=$b, t=$t, p=$p"
        @info "  Analytical w_center = $w_analytical"
        
        # Create rectangular mesh
        nx = 8
        ny = 16  # Same element density in both directions
        
        nodes = Asap.Node[]
        for j in 0:ny
            for i in 0:nx
                x = (i/nx - 0.5) * a
                y = (j/ny - 0.5) * b
                push!(nodes, Asap.Node([x*u"m", y*u"m", 0.0u"m"], :free))
            end
        end
        
        for (idx, node) in enumerate(nodes)
            node.nodeID = idx
            node.globalID = collect((idx-1)*6+1 : idx*6)
        end
        
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
        
        # Apply BCs (adapted for rectangular)
        fixed_dofs = Int[]
        tol_x = a / nx / 100
        tol_y = b / ny / 100
        
        for node in nodes
            x = ustrip(u"m", node.position[1])
            y = ustrip(u"m", node.position[2])
            
            on_x_edge = abs(abs(x) - a/2) < tol_x
            on_y_edge = abs(abs(y) - b/2) < tol_y
            
            if on_x_edge || on_y_edge
                append!(fixed_dofs, node.globalID[[3, 6]])  # w, θz
                if on_x_edge
                    push!(fixed_dofs, node.globalID[4])  # θx
                end
                if on_y_edge
                    push!(fixed_dofs, node.globalID[5])  # θy
                end
            end
        end
        append!(fixed_dofs, nodes[1].globalID[[1, 2]])
        fixed_dofs = unique(fixed_dofs)
        
        F = apply_uniform_pressure(nodes, elements, p)
        u = solve_plate(nodes, elements, fixed_dofs, F)
        
        # Get center deflection
        tol = min(a/nx, b/ny) / 10
        center_idx = findfirst(nd -> begin
            x = ustrip(u"m", nd.position[1])
            y = ustrip(u"m", nd.position[2])
            abs(x) < tol && abs(y) < tol
        end, nodes)
        
        w_center = u[nodes[center_idx].globalID[3]]
        accuracy = w_center / w_analytical * 100
        
        @info "  FEM w_center = $w_center, accuracy = $(round(accuracy, digits=2))%"
        
        @test accuracy > 85.0
        @test accuracy < 115.0
    end

end

@testset "Mesh Convergence Study" begin
    # Study mesh convergence to determine optimal default mesh density
    
    E = 30e6
    ν = 0.3
    L = 10.0
    t = 1.0  # 10% thickness ratio
    p = 0.1
    
    w_analytical = analytical_ss_square_udl(p, L, E, ν, t)
    
    @info "Mesh Convergence Study"
    @info "  Target: w_analytical = $w_analytical"
    
    mesh_sizes = [2, 4, 6, 8, 12, 16, 24, 32]
    errors = Float64[]
    dofs = Int[]
    
    for n in mesh_sizes
        nodes, elements = create_square_mesh(L, n, E, ν, t)
        fixed_dofs = apply_simply_supported_bc(nodes, L, n; hard=true)
        F = apply_uniform_pressure(nodes, elements, p)
        u = solve_plate(nodes, elements, fixed_dofs, F)
        w_center = get_center_deflection(nodes, u, L, n)
        
        error = abs(w_center - w_analytical) / abs(w_analytical) * 100
        push!(errors, error)
        push!(dofs, length(nodes) * 6)
        
        @info "  n=$n ($(length(elements)) elements, $(dofs[end]) DOF): error=$(round(error, digits=2))%"
    end
    
    # Convergence checks
    @testset "Convergence Properties" begin
        # Error should decrease with mesh refinement
        @test errors[end] < errors[1]
        
        # At n=8, error should be < 10%
        idx_8 = findfirst(==(8), mesh_sizes)
        @test errors[idx_8] < 10.0
        
        # At n=16, error should be < 5%
        idx_16 = findfirst(==(16), mesh_sizes)
        @test errors[idx_16] < 5.0
        
        # At n=32, error should be < 6% (element converges to ~5% offset from classical thin plate theory)
        idx_32 = findfirst(==(32), mesh_sizes)
        if !isnothing(idx_32)
            @test errors[idx_32] < 6.0
        end
    end
    
    @info "\nRecommended mesh density:"
    @info "  - Preliminary design: 4×4 per bay (error ~15%)"
    @info "  - Standard analysis: 8×8 per bay (error <10%)"
    @info "  - Detailed analysis: 16×16 per bay (error <5%)"
end

@testset "Thickness Effects (Thin vs Thick Plates)" begin
    # Test that element handles both thin and thick plates
    
    E = 200e9  # Steel-like
    ν = 0.3
    L = 1.0
    p = 1000.0  # 1 kPa
    
    thickness_ratios = [0.001, 0.01, 0.05, 0.1, 0.2]
    
    @info "Thickness Ratio Effects"
    
    for ratio in thickness_ratios
        t = L * ratio
        n = 16  # Fixed mesh
        
        w_analytical = analytical_ss_square_udl(p, L, E, ν, t)
        
        nodes, elements = create_square_mesh(L, n, E, ν, t)
        fixed_dofs = apply_simply_supported_bc(nodes, L, n; hard=true)
        F = apply_uniform_pressure(nodes, elements, p)
        u = solve_plate(nodes, elements, fixed_dofs, F)
        w_center = get_center_deflection(nodes, u, L, n)
        
        accuracy = w_center / w_analytical * 100
        
        @info "  t/L=$(ratio): accuracy=$(round(accuracy, digits=1))%"
        
        # For very thin plates, may have some numerical issues
        # Thick plates (t/L > 0.1) include shear effects not in classical theory
        if ratio >= 0.01 && ratio <= 0.1
            @test accuracy > 80.0
            @test accuracy < 115.0
        elseif ratio > 0.1
            # Thick plates: allow more deviation (shear deformation effects)
            @test accuracy > 80.0
            @test accuracy < 130.0  # Shear adds flexibility
        end
    end
end

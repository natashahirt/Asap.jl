#=
Graded/Structured Mesh Tests + Visual Inspection
=================================================

Tests for target_edge_length, refinement_edge_length, and auto-detected 
refinement targets. Generates PNG plots for visual inspection.
=#

using Test
using Asap
using Unitful

# ============================================================================
# Helper: extract mesh wireframe for plotting
# ============================================================================

"""Extract (x, y) arrays for every triangle edge — ready for plotting."""
function _extract_wireframe(shells)
    xs = Float64[]
    ys = Float64[]
    for s in shells
        nodes = s.nodes
        for (a, b) in ((1,2), (2,3), (3,1))
            push!(xs, ustrip(u"m", nodes[a].position[1]))
            push!(xs, ustrip(u"m", nodes[b].position[1]))
            push!(xs, NaN)
            push!(ys, ustrip(u"m", nodes[a].position[2]))
            push!(ys, ustrip(u"m", nodes[b].position[2]))
            push!(ys, NaN)
        end
    end
    return xs, ys
end

"""Extract unique node positions."""
function _extract_node_positions(shells)
    seen = Set{UInt64}()
    xs = Float64[]
    ys = Float64[]
    for s in shells, nd in s.nodes
        id = objectid(nd)
        id in seen && continue
        push!(seen, id)
        push!(xs, ustrip(u"m", nd.position[1]))
        push!(ys, ustrip(u"m", nd.position[2]))
    end
    return xs, ys
end

# ============================================================================
# Unit tests
# ============================================================================

@testset "Structured Mesh Generation" begin
    
    @testset "_t3block uniform" begin
        pts, conn = Asap._t3block(4.0, 3.0, 4, 3)
        @test length(pts) == 5 * 4  # (4+1)*(3+1)
        @test length(conn) == 2 * 4 * 3  # 2 tri per cell
        
        # Check bounds
        @test minimum(p[1] for p in pts) ≈ 0.0
        @test maximum(p[1] for p in pts) ≈ 4.0
        @test minimum(p[2] for p in pts) ≈ 0.0
        @test maximum(p[2] for p in pts) ≈ 3.0
    end
    
    @testset "_t3blockx graded" begin
        xs = [0.0, 0.1, 0.3, 0.6, 1.0, 2.0, 3.0, 4.0]
        ys = [0.0, 0.5, 1.0, 2.0, 3.0]
        pts, conn = Asap._t3blockx(xs, ys)
        @test length(pts) == length(xs) * length(ys)
        @test length(conn) == 2 * (length(xs)-1) * (length(ys)-1)
    end
    
    @testset "_graded_spacing" begin
        sp = Asap._graded_spacing(0.0, 10.0, [5.0], 1.0, 0.2, 2.0)
        
        # Should contain endpoints and target
        @test sp[1] ≈ 0.0
        @test sp[end] ≈ 10.0
        @test 5.0 in sp
        
        # Spacing near target should be small (≤ ~2× h_near)
        idx = findfirst(==(5.0), sp)
        @test idx !== nothing
        if idx !== nothing && idx > 1
            @test sp[idx] - sp[idx-1] ≤ 0.5  # graded, so last step before target ≤ ~2×h_near
        end
        
        # Multiple targets
        sp2 = Asap._graded_spacing(0.0, 10.0, [2.0, 8.0], 1.0, 0.15, 1.5)
        @test 2.0 in sp2
        @test 8.0 in sp2
    end
    
    @testset "_uniform_spacing" begin
        sp = Asap._uniform_spacing(0.0, 6.0, 1.5)
        @test sp[1] ≈ 0.0
        @test sp[end] ≈ 6.0
        @test length(sp) == 5  # ceil(6/1.5) + 1
    end
    
    @testset "_is_rectangular" begin
        rect = [(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]
        @test Asap._is_rectangular(rect) == true
        
        tri = [(0.0, 0.0), (4.0, 0.0), (2.0, 3.0)]
        @test Asap._is_rectangular(tri) == false
        
        non_rect = [(0.0, 0.0), (4.0, 0.0), (5.0, 3.0), (0.0, 3.0)]
        @test Asap._is_rectangular(non_rect) == false
    end
end

@testset "Shell with target_edge_length" begin
    sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
    
    @testset "Uniform structured mesh" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([6.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        shells = Shell((n1, n2, n3, n4), sec; target_edge_length=1.0u"m")
        @test length(shells) > 0
        
        # Structured mesh on 6×4 with h=1.0 → 6*4 cells * 2 tri = 48
        @test length(shells) == 48
        
        # Corner nodes preserved
        nodes = get_nodes(shells)
        @test n1 in nodes
        @test n2 in nodes
    end
    
    @testset "Structured finer" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        shells_coarse = Shell((n1, n2, n3, n4), sec; target_edge_length=2.0u"m")
        
        n1b = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2b = Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3b = Node([4.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        n4b = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        shells_fine = Shell((n1b, n2b, n3b, n4b), sec; target_edge_length=0.5u"m")
        
        @test length(shells_fine) > length(shells_coarse)
    end
    
    @testset "Non-rectangular falls back to Delaunay" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([3.0u"m", 3.0u"m", 0.0u"m"], :pinned)
        n4 = Node([1.0u"m", 4.0u"m", 0.0u"m"], :pinned)
        
        shells = Shell((n1, n2, n3, n4), sec; target_edge_length=0.5u"m")
        @test length(shells) > 0
        
        nodes = get_nodes(shells)
        @test n1 in nodes
    end
end

@testset "Shell with graded refinement" begin
    sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
    
    @testset "Rectangular + interior column" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        col = Node([3.0u"m", 3.0u"m", 0.0u"m"], :pinned)
        
        shells = Shell((n1, n2, n3, n4), sec;
            target_edge_length = 1.0u"m",
            interior_nodes = [col],
            refinement_edge_length = 0.2u"m",
            edge_support_type = :pinned
        )
        
        @test length(shells) > 48  # more than uniform 6×6 because of grading
        nodes = get_nodes(shells)
        @test col in nodes  # column node is structurally connected
    end
    
    @testset "Explicit refinement_targets" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        
        target_pt = Node([2.0u"m", 2.0u"m", 0.0u"m"], :free)
        
        shells = Shell((n1, n2, n3, n4), sec;
            target_edge_length = 1.0u"m",
            refinement_edge_length = 0.2u"m",
            refinement_targets = [target_pt]
        )
        @test length(shells) > 0
    end
    
    @testset "Delaunay + refinement rings" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([4.0u"m", 5.0u"m", 0.0u"m"], :pinned)
        n4 = Node([1.0u"m", 5.0u"m", 0.0u"m"], :pinned)
        col = Node([3.0u"m", 2.5u"m", 0.0u"m"], :pinned)
        
        shells_base = Shell((n1, n2, n3, n4), sec; target_edge_length=1.0u"m",
            interior_nodes=[col], edge_support_type=:free)
        shells_refined = Shell(
            (Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned),
             Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned),
             Node([4.0u"m", 5.0u"m", 0.0u"m"], :pinned),
             Node([1.0u"m", 5.0u"m", 0.0u"m"], :pinned)),
            sec;
            target_edge_length = 1.0u"m",
            interior_nodes = [Node([3.0u"m", 2.5u"m", 0.0u"m"], :pinned)],
            refinement_edge_length = 0.2u"m",
            edge_support_type = :free
        )
        
        @test length(shells_refined) > length(shells_base)
    end
end

@testset "Auto-detect refinement targets" begin
    sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
    
    @testset "Pinned corners become targets" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
        n3 = Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        n4 = Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
        
        # With refinement but no explicit targets → auto-detect from pinned corners
        shells = Shell((n1, n2, n3, n4), sec;
            target_edge_length = 1.0u"m",
            refinement_edge_length = 0.2u"m",
            edge_support_type = :pinned  # corners are supported
        )
        
        @test length(shells) > 48  # more than uniform because graded near corners
    end
    
    @testset "Free edges — no corner refinement" begin
        n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :free)
        n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :free)
        n3 = Node([6.0u"m", 6.0u"m", 0.0u"m"], :free)
        n4 = Node([0.0u"m", 6.0u"m", 0.0u"m"], :free)
        col = Node([3.0u"m", 3.0u"m", 0.0u"m"], :pinned)
        
        # Free edges + interior column → only refine around column
        shells = Shell((n1, n2, n3, n4), sec;
            target_edge_length = 1.0u"m",
            interior_nodes = [col],
            refinement_edge_length = 0.2u"m",
            edge_support_type = :free
        )
        @test length(shells) > 0
        @test col in get_nodes(shells)
    end
end

# ============================================================================
# Visual Inspection Plots (saved as PNG)
# ============================================================================

# Try loading a Makie backend for plotting; skip if unavailable
plot_ok = try
    @eval using CairoMakie
    true
catch
    try
        @eval using GLMakie
        true
    catch
        @warn "No Makie backend available — skipping visual mesh plots"
        false
    end
end

if plot_ok
    @testset "Visual mesh inspection" begin
        sec = ShellSection(0.15u"m", 30u"GPa", 0.2)
        output_dir = joinpath(@__DIR__, "..", "..", "..", "test_outputs")
        mkpath(output_dir)
        
        # ── 1. Uniform structured mesh ──
        let
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n3 = Node([6.0u"m", 4.0u"m", 0.0u"m"], :pinned)
            n4 = Node([0.0u"m", 4.0u"m", 0.0u"m"], :pinned)
            
            shells = Shell((n1, n2, n3, n4), sec; target_edge_length=0.5u"m")
            ex, ey = _extract_wireframe(shells)
            nx, ny = _extract_node_positions(shells)
            
            fig = Figure(size=(800, 600))
            ax = Axis(fig[1,1]; title="1) Uniform Structured (h=0.5m, 6×4)", 
                       xlabel="x [m]", ylabel="y [m]", aspect=DataAspect())
            lines!(ax, ex, ey; color=:steelblue, linewidth=0.5)
            scatter!(ax, nx, ny; color=:black, markersize=3)
            
            save(joinpath(output_dir, "mesh_01_uniform_structured.png"), fig)
            @test true
        end
        
        # ── 2. Graded structured mesh (1 interior column) ──
        let
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n3 = Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
            n4 = Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
            col = Node([3.0u"m", 3.0u"m", 0.0u"m"], :pinned)
            
            shells = Shell((n1, n2, n3, n4), sec;
                target_edge_length = 1.0u"m",
                interior_nodes = [col],
                refinement_edge_length = 0.15u"m",
                edge_support_type = :pinned
            )
            ex, ey = _extract_wireframe(shells)
            nx, ny = _extract_node_positions(shells)
            
            fig = Figure(size=(800, 800))
            ax = Axis(fig[1,1]; title="2) Delaunay + Rings — 1 column at (3,3), h_far=1.0, h_near=0.15",
                       xlabel="x [m]", ylabel="y [m]", aspect=DataAspect())
            lines!(ax, ex, ey; color=:steelblue, linewidth=0.5)
            scatter!(ax, nx, ny; color=:black, markersize=2)
            scatter!(ax, [3.0], [3.0]; color=:red, markersize=10, marker=:xcross)
            
            save(joinpath(output_dir, "mesh_02_graded_1col.png"), fig)
            @test true
        end
        
        # ── 3. Graded structured mesh (4 interior columns) ──
        let
            n1 = Node([0.0u"m",  0.0u"m",  0.0u"m"], :pinned)
            n2 = Node([12.0u"m", 0.0u"m",  0.0u"m"], :pinned)
            n3 = Node([12.0u"m", 12.0u"m", 0.0u"m"], :pinned)
            n4 = Node([0.0u"m",  12.0u"m", 0.0u"m"], :pinned)
            c1 = Node([4.0u"m",  4.0u"m",  0.0u"m"], :pinned)
            c2 = Node([8.0u"m",  4.0u"m",  0.0u"m"], :pinned)
            c3 = Node([8.0u"m",  8.0u"m",  0.0u"m"], :pinned)
            c4 = Node([4.0u"m",  8.0u"m",  0.0u"m"], :pinned)
            
            shells = Shell((n1, n2, n3, n4), sec;
                target_edge_length = 1.5u"m",
                interior_nodes = [c1, c2, c3, c4],
                refinement_edge_length = 0.2u"m",
                edge_support_type = :pinned
            )
            ex, ey = _extract_wireframe(shells)
            nx, ny = _extract_node_positions(shells)
            
            fig = Figure(size=(800, 800))
            ax = Axis(fig[1,1]; title="3) Delaunay + Rings — 4 columns, h_far=1.5, h_near=0.2",
                       xlabel="x [m]", ylabel="y [m]", aspect=DataAspect())
            lines!(ax, ex, ey; color=:steelblue, linewidth=0.4)
            scatter!(ax, nx, ny; color=:black, markersize=1.5)
            scatter!(ax, [4, 8, 8, 4], [4, 4, 8, 8]; color=:red, markersize=10, marker=:xcross)
            
            save(joinpath(output_dir, "mesh_03_graded_4col.png"), fig)
            @test true
        end
        
        # ── 4. Delaunay with refinement rings (irregular polygon) ──
        let
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n2 = Node([8.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n3 = Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
            n4 = Node([2.0u"m", 7.0u"m", 0.0u"m"], :pinned)
            col = Node([4.0u"m", 3.0u"m", 0.0u"m"], :pinned)
            
            shells = Shell((n1, n2, n3, n4), sec;
                target_edge_length = 1.0u"m",
                interior_nodes = [col],
                refinement_edge_length = 0.2u"m",
                edge_support_type = :free
            )
            ex, ey = _extract_wireframe(shells)
            nx, ny = _extract_node_positions(shells)
            
            fig = Figure(size=(800, 800))
            ax = Axis(fig[1,1]; title="4) Delaunay + Refinement Rings — irregular polygon",
                       xlabel="x [m]", ylabel="y [m]", aspect=DataAspect())
            lines!(ax, ex, ey; color=:steelblue, linewidth=0.5)
            scatter!(ax, nx, ny; color=:black, markersize=2)
            scatter!(ax, [4.0], [3.0]; color=:red, markersize=10, marker=:xcross)
            
            save(joinpath(output_dir, "mesh_04_delaunay_refined.png"), fig)
            @test true
        end
        
        # ── 5. Legacy n-based mesh (backward compatibility) ──
        let
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n2 = Node([4.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n3 = Node([4.0u"m", 3.0u"m", 0.0u"m"], :pinned)
            n4 = Node([0.0u"m", 3.0u"m", 0.0u"m"], :pinned)
            
            shells = Shell((n1, n2, n3, n4), 6, sec)
            ex, ey = _extract_wireframe(shells)
            nx, ny = _extract_node_positions(shells)
            
            fig = Figure(size=(800, 600))
            ax = Axis(fig[1,1]; title="5) Legacy n=6 (Delaunay, backward compatible)",
                       xlabel="x [m]", ylabel="y [m]", aspect=DataAspect())
            lines!(ax, ex, ey; color=:steelblue, linewidth=0.5)
            scatter!(ax, nx, ny; color=:black, markersize=3)
            
            save(joinpath(output_dir, "mesh_05_legacy_n6.png"), fig)
            @test true
        end
        
        # ── 6. Comparison: with vs without refinement ──
        let
            fig = Figure(size=(1400, 600))
            
            # Without refinement
            n1 = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n2 = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n3 = Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
            n4 = Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
            col1 = Node([3.0u"m", 3.0u"m", 0.0u"m"], :pinned)
            
            shells_no_ref = Shell((n1, n2, n3, n4), sec;
                target_edge_length = 0.5u"m",
                interior_nodes = [col1],
                edge_support_type = :pinned
            )
            
            ax1 = Axis(fig[1,1]; title="Without refinement (h=0.5m uniform)",
                        xlabel="x [m]", ylabel="y [m]", aspect=DataAspect())
            ex, ey = _extract_wireframe(shells_no_ref)
            lines!(ax1, ex, ey; color=:steelblue, linewidth=0.5)
            scatter!(ax1, [3.0], [3.0]; color=:red, markersize=10, marker=:xcross)
            
            # With refinement
            n1b = Node([0.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n2b = Node([6.0u"m", 0.0u"m", 0.0u"m"], :pinned)
            n3b = Node([6.0u"m", 6.0u"m", 0.0u"m"], :pinned)
            n4b = Node([0.0u"m", 6.0u"m", 0.0u"m"], :pinned)
            col2 = Node([3.0u"m", 3.0u"m", 0.0u"m"], :pinned)
            
            shells_ref = Shell((n1b, n2b, n3b, n4b), sec;
                target_edge_length = 0.5u"m",
                interior_nodes = [col2],
                refinement_edge_length = 0.1u"m",
                edge_support_type = :pinned
            )
            
            ax2 = Axis(fig[1,2]; title="With refinement — Delaunay (h_far=0.5, h_near=0.1)",
                        xlabel="x [m]", ylabel="y [m]", aspect=DataAspect())
            ex2, ey2 = _extract_wireframe(shells_ref)
            lines!(ax2, ex2, ey2; color=:steelblue, linewidth=0.5)
            scatter!(ax2, [3.0], [3.0]; color=:red, markersize=10, marker=:xcross)
            
            save(joinpath(output_dir, "mesh_06_comparison.png"), fig)
            @test length(shells_ref) > length(shells_no_ref)
        end
        
        # ── 7. Random polygon with random interior columns ──
        let
            using Random
            rng = MersenneTwister(42)  # fixed seed for reproducibility
            
            # Generate a random convex polygon (sample angles, sort, then place on unit circle and scale)
            n_corners = rand(rng, 5:8)
            angles = sort!(rand(rng, n_corners) .* 2π)
            scale = 5.0 + rand(rng) * 10.0  # 5–15 m
            cx, cy = 7.0, 7.0  # center
            corner_xs = cx .+ scale .* cos.(angles)
            corner_ys = cy .+ scale .* sin.(angles)
            
            corners = tuple([Node([corner_xs[i]*u"m", corner_ys[i]*u"m", 0.0u"m"], :pinned) 
                             for i in 1:n_corners]...)
            
            # Random interior columns (ensure they're inside the convex hull)
            n_cols = rand(rng, 2:5)
            col_nodes = Node[]
            for _ in 1:n_cols
                # Random convex combination of corners → guaranteed interior
                weights = rand(rng, n_corners)
                weights ./= sum(weights)
                px = sum(weights .* corner_xs)
                py = sum(weights .* corner_ys)
                push!(col_nodes, Node([px*u"m", py*u"m", 0.0u"m"], :pinned))
            end
            
            shells = Shell(corners, sec;
                target_edge_length = 1.0u"m",
                interior_nodes = col_nodes,
                refinement_edge_length = 0.2u"m",
                edge_support_type = :free
            )
            
            ex, ey = _extract_wireframe(shells)
            nx, ny = _extract_node_positions(shells)
            col_x = [ustrip(u"m", c.position[1]) for c in col_nodes]
            col_y = [ustrip(u"m", c.position[2]) for c in col_nodes]
            
            fig = Figure(size=(800, 800))
            ax = Axis(fig[1,1]; 
                title="7) Random $(n_corners)-gon, $(n_cols) random columns (seed=42)",
                xlabel="x [m]", ylabel="y [m]", aspect=DataAspect())
            lines!(ax, ex, ey; color=:steelblue, linewidth=0.5)
            scatter!(ax, nx, ny; color=:black, markersize=2)
            scatter!(ax, col_x, col_y; color=:red, markersize=10, marker=:xcross)
            # Draw polygon outline
            poly_x = [corner_xs; corner_xs[1]]
            poly_y = [corner_ys; corner_ys[1]]
            lines!(ax, poly_x, poly_y; color=:red, linewidth=1.5, linestyle=:dash)
            
            save(joinpath(output_dir, "mesh_07_random_polygon.png"), fig)
            @test length(shells) > 0
            # All column nodes should be in the mesh
            mesh_nodes = get_nodes(shells)
            for c in col_nodes
                @test c in mesh_nodes
            end
        end
        
        @info "Mesh plots saved to: $output_dir"
    end
end

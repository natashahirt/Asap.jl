# =============================================================================
# Visual Debug: Cell Depth Computation for EFM l2
# Run: julia --project=. external/Asap/test/tributary/debug_cell_depths.jl
# =============================================================================

using Pkg
Pkg.activate(".")

using Revise
using Asap
using Meshes
using GLMakie
using Unitful: @u_str

# =============================================================================
# Test Shape Definitions
# =============================================================================

"""Create Meshes.Point vertices from (x,y) tuples (in meters)."""
function make_vertices(coords::Vector{NTuple{2,Float64}})
    return [Meshes.Point(c[1] * u"m", c[2] * u"m") for c in coords]
end

# --- Regular Shapes ---
rectangle() = make_vertices([(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0)])
square() = make_vertices([(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])
parallelogram() = make_vertices([(0.0, 0.0), (5.0, 0.0), (6.5, 3.0), (1.5, 3.0)])
trapezoid() = make_vertices([(0.0, 0.0), (6.0, 0.0), (5.0, 3.0), (1.0, 3.0)])
trapezoid_wide_top() = make_vertices([(1.0, 0.0), (5.0, 0.0), (7.0, 4.0), (-1.0, 4.0)])
pentagon() = make_vertices([(2.0, 0.0), (4.0, 0.0), (5.0, 2.5), (3.0, 4.0), (1.0, 2.5)])
hexagon() = make_vertices([(1.0, 0.0), (3.0, 0.0), (4.0, 1.5), (3.0, 3.0), (1.0, 3.0), (0.0, 1.5)])

# --- Irregular / Concave Shapes ---
irregular_quad() = make_vertices([(0.0, 0.0), (5.0, 0.5), (4.5, 3.5), (0.5, 2.5)])
l_shape() = make_vertices([(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 4.0), (0.0, 4.0)])
arrow_shape() = make_vertices([(2.0, 0.0), (4.0, 2.0), (3.0, 2.0), (3.0, 4.0), (1.0, 4.0), (1.0, 2.0), (0.0, 2.0)])
chevron() = make_vertices([(0.0, 0.0), (2.0, 2.0), (4.0, 0.0), (4.0, 1.0), (2.0, 3.0), (0.0, 1.0)])
long_thin_rect() = make_vertices([(0.0, 0.0), (10.0, 0.0), (10.0, 2.0), (0.0, 2.0)])
narrow_triangle() = make_vertices([(0.0, 0.0), (8.0, 0.0), (4.0, 2.0)])

# =============================================================================
# Visualization
# =============================================================================

"""Convert Meshes.Point to (x,y) tuple."""
_to_2d(v::Meshes.Point) = (Float64(Meshes.coords(v).x.val), Float64(Meshes.coords(v).y.val))

"""Shoelace formula for signed polygon area."""
function _polygon_area(pts::Vector{NTuple{2,Float64}})
    n = length(pts)
    n < 3 && return 0.0
    sum(pts[i][1] * pts[mod1(i+1,n)][2] - pts[mod1(i+1,n)][1] * pts[i][2] for i in 1:n) / 2
end

"""Ensure polygon is CCW oriented."""
_ensure_ccw(pts::Vector{NTuple{2,Float64}}) = _polygon_area(pts) > 0 ? pts : reverse(pts)

"""Check if point is inside polygon (ray casting)."""
function _point_in_polygon(pt, poly)
    x, y = pt
    n = length(poly)
    inside = false
    j = n
    for i in 1:n
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

"""Plot a shape with all sampled depth vectors and l2_stiff."""
function plot_with_depths!(ax, verts; title="", show_tribs=true)
    colors = [:coral, :skyblue, :lightgreen, :plum, :gold, :salmon, :cyan, :pink]
    
    # Use same canonical CCW ordering as computation
    pts = [_to_2d(v) for v in verts]
    pts = _ensure_ccw(pts)
    n = length(pts)
    
    # Center for display
    cx = sum(p[1] for p in pts) / n
    cy = sum(p[2] for p in pts) / n
    pts_centered = [(p[1] - cx, p[2] - cy) for p in pts]
    
    results = get_tributary_polygons_isotropic(verts)
    
    # Layer 1: Tributary polygons (background)
    if show_tribs
        for (i, trib) in enumerate(results)
            if !isempty(trib.s)
                beam_start = pts[trib.local_edge_idx]
                beam_end = pts[mod1(trib.local_edge_idx + 1, n)]
                abs_verts = Asap.vertices(trib, beam_start, beam_end)
                
                if !isempty(abs_verts)
                    txs = [v[1] - cx for v in abs_verts]
                    tys = [v[2] - cy for v in abs_verts]
                    push!(txs, txs[1]); push!(tys, tys[1])
                    
                    c = colors[mod1(i, length(colors))]
                    poly!(ax, Point2f.(txs, tys), color=(c, 0.2), strokecolor=(c, 0.4), strokewidth=0.5)
                end
            end
        end
    end
    
    # Layer 2: Shape vertices
    scatter!(ax, Point2f.(pts_centered), color=:black, markersize=6)
    
    # Layer 3: Cell depth vectors (all samples + max highlighted)
    for trib in results
        edge_idx = trib.local_edge_idx
        v1 = pts_centered[edge_idx]
        v2 = pts_centered[mod1(edge_idx + 1, n)]
        
        edge_vec = (v2[1] - v1[1], v2[2] - v1[2])
        edge_len = hypot(edge_vec...)
        edge_len < 1e-9 && continue
        
        edge_dir = (edge_vec[1] / edge_len, edge_vec[2] / edge_len)
        perp_dir = (-edge_dir[2], edge_dir[1])
        
        point_at_s(s) = (v1[1] + s * (v2[1] - v1[1]), v1[2] + s * (v2[2] - v1[2]))
        mid = point_at_s(0.5)
        
        # Ensure perp_dir points inward
        test_pt = (mid[1] + 1e-4 * perp_dir[1], mid[2] + 1e-4 * perp_dir[2])
        if !_point_in_polygon(test_pt, pts_centered)
            perp_dir = (-perp_dir[1], -perp_dir[2])
        end
        
        # Draw ALL sampled depth vectors (gray, thin)
        for (s_val, depth) in zip(trib.cell_depths_s, trib.cell_depths)
            pt = point_at_s(s_val)
            end_pt = (pt[1] + perp_dir[1] * depth, pt[2] + perp_dir[2] * depth)
            lines!(ax, [pt[1], end_pt[1]], [pt[2], end_pt[2]], 
                   color=(:gray, 0.5), linewidth=1)
            scatter!(ax, Point2f(pt), color=:gray, markersize=3)
        end
        
        # Max depth - red dashed (highlighted)
        pt_max = point_at_s(trib.cell_depth_s_max)
        end_max = (pt_max[1] + perp_dir[1] * trib.cell_depth_max, 
                   pt_max[2] + perp_dir[2] * trib.cell_depth_max)
        lines!(ax, [pt_max[1], end_max[1]], [pt_max[2], end_max[2]], 
               color=(:red, 0.8), linewidth=2, linestyle=:dash)
        scatter!(ax, Point2f(pt_max), color=:red, markersize=5, marker=:diamond)
        
        # Label with l2_stiff (cubic mean)
        label_pos = (mid[1] + perp_dir[1] * (trib.l2_stiff * 0.5 + 0.3), 
                     mid[2] + perp_dir[2] * (trib.l2_stiff * 0.5 + 0.3))
        text!(ax, label_pos[1], label_pos[2], text="$(round(trib.l2_stiff, digits=2))m", 
              fontsize=8, align=(:center, :center), color=:blue)
    end
    
    ax.title = title
end

"""Create the full visualization."""
function visualize_cell_depths()
    shapes = [
        ("Square (4×4)", square()),
        ("Rectangle (6×4)", rectangle()),
        ("Long Rectangle (10×2)", long_thin_rect()),
        ("Trapezoid", trapezoid()),
        ("Trapezoid (wide top)", trapezoid_wide_top()),
        ("Parallelogram", parallelogram()),
        ("Pentagon", pentagon()),
        ("Hexagon", hexagon()),
        ("Triangle", narrow_triangle()),
        ("Irregular Quad", irregular_quad()),
        ("L-Shape", l_shape()),
        ("Arrow", arrow_shape()),
        ("Chevron", chevron()),
    ]
    
    n_shapes = length(shapes)
    n_cols = 4
    n_rows = ceil(Int, n_shapes / n_cols)
    
    fig = Figure(size = (400 * n_cols, 400 * n_rows), fontsize=12)
    
    for (i, (name, verts)) in enumerate(shapes)
        row = div(i - 1, n_cols) + 1
        col = mod(i - 1, n_cols) + 1
        ax = Axis(fig[row, col], aspect=DataAspect(), xlabel="x [m]", ylabel="y [m]")
        plot_with_depths!(ax, verts; title=name)
    end
    
    # Legend
    legend_ax = Axis(fig[n_rows + 1, 1:n_cols], height=50)
    hidedecorations!(legend_ax); hidespines!(legend_ax)
    text!(legend_ax, 0.1, 0.5, text="Legend: ", fontsize=14, font=:bold)
    lines!(legend_ax, [0.25, 0.35], [0.5, 0.5], color=:gray, linewidth=1)
    text!(legend_ax, 0.36, 0.5, text="sampled depths", fontsize=12, color=:gray)
    lines!(legend_ax, [0.55, 0.65], [0.5, 0.5], color=:red, linewidth=2, linestyle=:dash)
    text!(legend_ax, 0.66, 0.5, text="max depth", fontsize=12, color=:red)
    text!(legend_ax, 0.85, 0.5, text="labels = l2_stiff (cubic mean)", fontsize=12, color=:blue)
    
    Label(fig[0, :], "Cell Depths & l2_stiff (Stiffness-Consistent Width for EFM)", fontsize=20, font=:bold)
    
    return fig
end

"""Print depth values for verification."""
function print_depth_values()
    shapes = [
        ("Square (4×4)", square()),
        ("Rectangle (6×4)", rectangle()),
        ("Trapezoid", trapezoid()),
        ("L-Shape", l_shape()),
        ("Irregular Quad", irregular_quad()),
        ("Chevron", chevron()),
    ]
    
    println("=" ^ 70)
    println("Cell Depth Values (for EFM l2 calculation)")
    println("=" ^ 70)
    
    for (name, verts) in shapes
        pts = [_to_2d(v) for v in verts]
        pts = _ensure_ccw(pts)
        n_pts = length(pts)
        results = get_tributary_polygons_isotropic(verts)
        
        println("\n$(name) ($(length(verts)) vertices)")
        println("-" ^ 50)
        
        for trib in results
            edge_idx = trib.local_edge_idx
            v1 = pts[edge_idx]
            v2 = pts[mod1(edge_idx + 1, n_pts)]
            
            # Show sampled depths
            samples = join(["$(round(d, digits=2))@$(round(s, digits=2))" 
                           for (s, d) in zip(trib.cell_depths_s, trib.cell_depths)], ", ")
            
            println("  Edge $(edge_idx): ($(round(v1[1], digits=1)), $(round(v1[2], digits=1))) → ($(round(v2[1], digits=1)), $(round(v2[2], digits=1)))")
            println("    l2_stiff=$(round(trib.l2_stiff, digits=3))m  max=$(round(trib.cell_depth_max, digits=2))m")
            println("    samples: $samples")
        end
    end
    println("\n" * "=" ^ 70)
end

# =============================================================================
# Main
# =============================================================================

print_depth_values()

println("\nGenerating visualization...")
fig = visualize_cell_depths()
display(fig)

println("\n✓ Complete! Labels show l2_stiff (cubic mean effective width)")

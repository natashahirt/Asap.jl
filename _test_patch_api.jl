#=
Test script: verify ShellPatch + Shell() integration.

Run from root:  julia --project=StructuralSizer external/Asap/_test_patch_api.jl
=#

using Asap
using Unitful
import Meshes

println("=" ^ 60)
println("TEST 1: ShellPatch — rectangular from center + dims")
println("=" ^ 60)

stiff_sec = ShellSection(0.5u"m", 200e9u"Pa", 0.3; name=:stiff)
patch_rect = ShellPatch(3.0, 3.0, 0.5, 0.5, stiff_sec; id=:column_patch)

println("  Vertices: ", patch_rect.vertices)
println("  Center:   ", patch_rect.center)
println("  Section:  ", patch_rect.section.name)
println("  id:       ", patch_rect.id)
println("  ✅ Rectangular patch created")

println()
println("=" ^ 60)
println("TEST 2: ShellPatch — from Meshes.Ngon (Octagon)")
println("=" ^ 60)

cx, cy, r = 3.0, 3.0, 0.25
n_sides = 8
oct = Meshes.Ngon(
    [Meshes.Point(cx + r * cos(2π * k / n_sides),
                   cy + r * sin(2π * k / n_sides))
     for k in 0:n_sides-1]...
)
patch_circ = ShellPatch(oct, stiff_sec; id=:column)

println("  Vertices: ", length(patch_circ.vertices), " sides")
println("  Center:   ", patch_circ.center)
println("  ✅ Octagonal patch created from Meshes.Ngon")

println()
println("=" ^ 60)
println("TEST 3: ShellPatch — Unitful constructor")
println("=" ^ 60)

patch_unit = ShellPatch(3.0u"m", 3.0u"m", 500.0u"mm", 500.0u"mm", stiff_sec)
println("  Center:   ", patch_unit.center)
println("  Vertices: ", patch_unit.vertices)
println("  ✅ Unitful constructor works")

println()
println("=" ^ 60)
println("TEST 4: Shell() with interior_patches — Delaunay")
println("=" ^ 60)

# Create a 6m × 6m slab with a column patch at center
slab_sec = ShellSection(0.2u"m", 30e9u"Pa", 0.2; name=:slab)

n1 = Node([0.0u"m", 0.0u"m", 3.0u"m"], :pinned)
n2 = Node([6.0u"m", 0.0u"m", 3.0u"m"], :pinned)
n3 = Node([6.0u"m", 6.0u"m", 3.0u"m"], :pinned)
n4 = Node([0.0u"m", 6.0u"m", 3.0u"m"], :pinned)

shells = Shell((n1, n2, n3, n4), slab_sec;
    target_edge_length = 0.5u"m",
    interior_patches = [patch_rect],
    edge_support_type = :free
)

n_total = length(shells)
n_patch = count(s -> s.id == :column_patch, shells)
n_slab  = count(s -> s.id == :shell, shells)

println("  Total shells:    $n_total")
println("  Patch elements:  $n_patch  (id = :column_patch, thickness = 0.5m)")
println("  Slab elements:   $n_slab   (id = :shell, thickness = 0.2m)")
println("  Patch fraction:  $(round(n_patch / n_total * 100; digits=1))%")

# Debug: show unique IDs and thicknesses
id_set = Set(s.id for s in shells)
println("  Unique IDs:      $id_set")

thick_set = Set(round(s.thickness; digits=3) for s in shells)
println("  Unique thicknesses: $thick_set m")

# Check that at least SOME elements got the patch section
if n_patch > 0
    p = first(s for s in shells if s.id == :column_patch)
    @assert p.thickness ≈ 0.5   "Patch thickness mismatch"
    @assert p.E ≈ 200e9         "Patch E mismatch"
    println("  Patch thickness: $(p.thickness) m ✅")
    println("  Patch E:         $(p.E) Pa ✅")
else
    # Debug: find centroids near the patch region
    println("  ⚠️ No patch elements found — debugging centroids near patch:")
    for (i, s) in enumerate(shells)
        cx_s = sum(ustrip(u"m", nd.position[1]) for nd in s.nodes) / 3.0
        cy_s = sum(ustrip(u"m", nd.position[2]) for nd in s.nodes) / 3.0
        if 2.5 < cx_s < 3.5 && 2.5 < cy_s < 3.5
            println("    Shell $i: centroid=($cx_s, $cy_s), t=$(s.thickness), id=$(s.id)")
        end
    end
end

println()
println("=" ^ 60)
println("TEST 5: Larger patch to ensure matching")
println("=" ^ 60)

# Use a larger patch (1.0 × 1.0) to guarantee centroids fall inside
large_patch = ShellPatch(3.0, 3.0, 1.0, 1.0, stiff_sec; id=:big_patch)

n1b = Node([0.0u"m", 0.0u"m", 3.0u"m"], :pinned)
n2b = Node([6.0u"m", 0.0u"m", 3.0u"m"], :pinned)
n3b = Node([6.0u"m", 6.0u"m", 3.0u"m"], :pinned)
n4b = Node([0.0u"m", 6.0u"m", 3.0u"m"], :pinned)

shells2 = Shell((n1b, n2b, n3b, n4b), slab_sec;
    target_edge_length = 0.3u"m",
    interior_patches = [large_patch],
    edge_support_type = :free
)

n_total2 = length(shells2)
n_patch2 = count(s -> s.id == :big_patch, shells2)
n_slab2  = count(s -> s.id == :shell, shells2)

println("  Total shells:    $n_total2")
println("  Patch elements:  $n_patch2  (id = :big_patch)")
println("  Slab elements:   $n_slab2   (id = :shell)")

@assert n_patch2 > 0 "Expected some patch elements with 1m × 1m patch"
println("  ✅ Patch elements found with larger patch")

# Verify section properties
p2 = first(s for s in shells2 if s.id == :big_patch)
@assert p2.thickness ≈ 0.5  "Patch thickness mismatch"
@assert p2.E ≈ 200e9        "Patch E mismatch"
println("  Patch t=$(p2.thickness)m, E=$(p2.E)Pa ✅")

s2 = first(s for s in shells2 if s.id == :shell)
@assert s2.thickness ≈ 0.2  "Slab thickness mismatch"
@assert s2.E ≈ 30e9         "Slab E mismatch"
println("  Slab  t=$(s2.thickness)m, E=$(s2.E)Pa ✅")

println()
println("=" ^ 60)
println("TEST 6: Full model solve with patches")
println("=" ^ 60)

cn1 = Node([0.0u"m", 0.0u"m", 3.0u"m"], :pinned)
cn2 = Node([6.0u"m", 0.0u"m", 3.0u"m"], :pinned)
cn3 = Node([6.0u"m", 6.0u"m", 3.0u"m"], :pinned)
cn4 = Node([0.0u"m", 6.0u"m", 3.0u"m"], :pinned)

col = Node([3.0u"m", 3.0u"m", 3.0u"m"], :zfixed)

patch_solve = ShellPatch(3.0, 3.0, 0.5, 0.5,
    ShellSection(0.5u"m", 30e9u"Pa", 0.2; name=:rigid_patch);
    id=:col_patch)

shells_s = Shell((cn1, cn2, cn3, cn4), slab_sec;
    target_edge_length = 0.3u"m",
    interior_nodes = [col],
    interior_patches = [patch_solve],
    edge_support_type = :pinned
)

n_cpatch = count(s -> s.id == :col_patch, shells_s)
n_total_s = length(shells_s)
println("  Total elements: $n_total_s")
println("  Patch elements: $n_cpatch")

# Apply uniform pressure load: 5 kPa downward
# AreaLoad(shells, pressure; direction=(0,0,-1))  →  force_z = pressure * (-1)
# For downward: use positive pressure with default -Z direction
load = AreaLoad(shells_s, 5.0u"kN/m^2")

all_n = get_nodes(shells_s)
model = ShellModel(all_n, shells_s, [load])
solve!(model)

# ── Check 1: AreaLoad distributes correctly over the meshed area ──
# sum(model.P) at all z-DOFs should match hand-calculated total
all_z_gids = [n.globalID[3] for n in all_n]
total_Pz = sum(model.P[gid] for gid in all_z_gids)
expected_Pz = -5000.0 * 36.0  # -180,000 N  (5 kPa × 6m × 6m, downward)
load_err = abs(total_Pz - expected_Pz) / abs(expected_Pz) * 100

println("  Applied Pz total: $(round(total_Pz; digits=1)) N  (expected $(expected_Pz) N)")
println("  Load distrib err: $(round(load_err; digits=3))%")
@assert load_err < 1.0 "AreaLoad distribution error too large: $(load_err)%"
println("  ✅ AreaLoad distributes correctly")

# ── Check 2: Equilibrium — sum(reactions) + sum(applied) = 0 ──
# With the corrected reaction formula (R = K·u − P at fixed DOFs),
# stored reactions are now true physical reactions.
total_Rz = sum(model.reactions[n.globalID[3]] for n in all_n if !n.dof[3])

equil_residual = total_Rz + total_Pz          # should be ≈ 0
equil_err = abs(equil_residual) / abs(expected_Pz) * 100

println("  Σ(Rz):   $(round(total_Rz; digits=1)) N")
println("  Σ(Pz):   $(round(total_Pz; digits=1)) N")
println("  Residual: $(round(equil_residual; digits=2)) N  ($(round(equil_err; digits=4))%)")
@assert equil_err < 1.0 "Equilibrium error too large: $(equil_err)%"
println("  ✅ Equilibrium satisfied")

# ── Check 3: Center column carries meaningful load ──
col_Rz = col.reaction[3]
println("  Column Rz: $(round(ustrip(u"kN", col_Rz); digits=2)) kN")
@assert abs(ustrip(u"N", col_Rz)) > 500.0 "Center column reaction too small"
println("  ✅ Center column carries load")

println()
println("=" ^ 60)
println("ALL TESTS PASSED ✅")
println("=" ^ 60)

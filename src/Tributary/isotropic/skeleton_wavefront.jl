"""Compute bisectors and speeds for active vertices (isotropic, all weights = 1)."""
function _compute_bisectors(pts::Vector{NTuple{2,Float64}}, n::Int)
    bisectors = Vector{NTuple{2,Float64}}(undef, n)
    speeds = Vector{Float64}(undef, n)
    
    for i in 1:n
        prev_i = mod1(i - 1, n)
        next_i = mod1(i + 1, n)
        
        p_prev = pts[prev_i]
        p_curr = pts[i]
        p_next = pts[next_i]
        
        v_in = (p_curr[1] - p_prev[1], p_curr[2] - p_prev[2])
        v_out = (p_next[1] - p_curr[1], p_next[2] - p_curr[2])
        
        len_in = hypot(v_in...)
        len_out = hypot(v_out...)
        
        if len_in < 1e-10 || len_out < 1e-10
            bisectors[i] = (0.0, 0.0)
            speeds[i] = 0.0
            continue
        end
        
        v_in = v_in ./ len_in
        v_out = v_out ./ len_out
        
        # Inward normals (90° CCW rotation)
        n_in = (-v_in[2], v_in[1])
        n_out = (-v_out[2], v_out[1])
        
        bx, by = n_in[1] + n_out[1], n_in[2] + n_out[2]
        blen = hypot(bx, by)
        
        if blen < 1e-10
            bisectors[i] = (0.0, 0.0)
            speeds[i] = 0.0
        else
            bisectors[i] = (bx / blen, by / blen)
            speeds[i] = 2.0 / blen
        end
    end
    
    return bisectors, speeds
end

# =============================================================================
# Isotropic wavefront with split events
# =============================================================================

"""
    _split_events(pts, bisectors, speeds) -> (t_min, events)

Split events of one CCW wavefront: a reflex vertex `i` reaching the moving
offset of a non-adjacent edge `j` (vertex `j` to `j + 1`) at time `t`.
`events` holds every `(i, j, t, hit)` with `t` at `t_min`, earliest first.
"""
function _split_events(pts::Vector{NTuple{2,Float64}}, bisectors, speeds)
    n = length(pts)
    found = Tuple{Int, Int, Float64, NTuple{2,Float64}}[]
    for i in 1:n
        a = pts[mod1(i - 1, n)]
        p = pts[i]
        c = pts[mod1(i + 1, n)]
        _cross(a, p, c) < 0 || continue
        speeds[i] > 0 || continue
        w = bisectors[i] .* speeds[i]
        for j in 1:n
            (j == i || j == mod1(i - 1, n)) && continue
            k = mod1(j + 1, n)
            ea, eb = pts[j], pts[k]
            L = _dist(ea, eb)
            L < 1e-12 && continue
            u = ((eb[1] - ea[1]) / L, (eb[2] - ea[2]) / L)
            nrm = (-u[2], u[1])
            d0 = (p[1] - ea[1]) * nrm[1] + (p[2] - ea[2]) * nrm[2]
            den = 1.0 - (w[1] * nrm[1] + w[2] * nrm[2])
            (d0 >= -1e-12 && den > 1e-12) || continue
            t = max(d0 / den, 0.0)
            hit = (p[1] + w[1] * t, p[2] + w[2] * t)
            # The offset edge at time t runs between its moving end vertices.
            at = pts[j] .+ (bisectors[j] .* (speeds[j] * t))
            bt = pts[k] .+ (bisectors[k] .* (speeds[k] * t))
            span = (bt[1] - at[1]) * u[1] + (bt[2] - at[2]) * u[2]
            s = (hit[1] - at[1]) * u[1] + (hit[2] - at[2]) * u[2]
            tol = 1e-9 * max(L, 1.0)
            (span > tol && -tol <= s <= span + tol) || continue
            push!(found, (i, j, t, hit))
        end
    end
    isempty(found) && return (Inf, found)
    t_min = minimum(e[3] for e in found)
    tol = max(1e-9, 1e-6 * t_min)
    events = sort!([e for e in found if e[3] <= t_min + tol]; by = e -> e[3])
    return (t_min, events)
end

"""
Drop a vertex where the boundary reverses along the same line.

The walk goes out to that vertex and comes straight back, so the vertex
encloses no area. The incoming edge keeps its face and runs on to the
vertex after the tip.
"""
function _drop_reversals(pts::Vector{NTuple{2,Float64}}, em::Vector{Int})
    pts = copy(pts)
    em = copy(em)
    changed = true
    while changed && length(pts) >= 4
        changed = false
        n = length(pts)
        for i in 1:n
            prev = pts[mod1(i - 1, n)]
            curr = pts[i]
            nxt = pts[mod1(i + 1, n)]
            vin = (curr[1] - prev[1], curr[2] - prev[2])
            vout = (nxt[1] - curr[1], nxt[2] - curr[2])
            lin, lout = hypot(vin...), hypot(vout...)
            (lin < 1e-9 || lout < 1e-9) && continue
            vin = vin ./ lin
            vout = vout ./ lout
            dotv = vin[1] * vout[1] + vin[2] * vout[2]
            # A tip that turns back along its incoming edge. A millimetre
            # kink still reverses; a real corner does not.
            dotv < -0.99 || continue
            chord = (nxt[1] - prev[1], nxt[2] - prev[2])
            Lc = hypot(chord...)
            dist = Lc < 1e-9 ? 0.0 :
                abs(chord[1] * (curr[2] - prev[2]) - chord[2] * (curr[1] - prev[1])) / Lc
            dist < 1e-3 || continue
            deleteat!(pts, i)
            deleteat!(em, i)
            changed = true
            break
        end
    end
    return pts, em
end

"""
Drop the zero-length edges of a wavefront. Removing edge `k` merges vertices
`k` and `k + 1`; the merged vertex sits between edges `k - 1` and `k + 1`.
"""
function _drop_tiny_edges(pts::Vector{NTuple{2,Float64}}, em::Vector{Int}; eps_len = 1e-9)
    pts = copy(pts)
    em = copy(em)
    k = 1
    while length(pts) >= 2 && k <= length(pts)
        if _dist(pts[k], pts[mod1(k + 1, length(pts))]) < eps_len
            # Vertex k coincides with k + 1, so it goes with edge k.
            deleteat!(em, k)
            deleteat!(pts, k)
            k = 1
        else
            k += 1
        end
    end
    return pts, em
end

"""
Record a two-vertex wavefront as the ridge between its two edges. Front edge
`em[1]` runs from `pts[1]` to `pts[2]` with its swept face on its right.
"""
function _record_ridge!(dcel, registry, pts, em)
    length(pts) == 2 || return
    em[1] == em[2] && return
    _dist(pts[1], pts[2]) > 1e-10 || return
    _record_skeleton_arc!(dcel, registry, pts[1], pts[2], em[2], em[1])
end

"""
    _propagate_isotropic!(dcel, registry, pts, edge_map; max_events)

Shrink the CCW polygon `pts` at unit speed and record the skeleton arcs in
`dcel`. `edge_map[k]` is the face of edge `k` (vertex `k` to `k + 1`).

Edge events remove an edge whose ends meet. Split events cut a front in two
where a reflex vertex reaches another edge of the same front; that edge then
lies in both pieces. Both kinds make the result the straight skeleton of a
concave polygon, not only of a convex one.
"""
function _propagate_isotropic!(dcel::DCEL, registry::VertexRegistry,
                               pts::Vector{NTuple{2,Float64}}, edge_map::Vector{Int};
                               max_events::Int = 50 * length(pts) + 100)
    pts, edge_map = _drop_reversals(pts, edge_map)
    length(pts) < 3 && return dcel
    # A stitched ring can come back clockwise. The shrink walks with the
    # slab on the left, so a clockwise front is reversed, edge faces with it.
    if _polygon_area(pts) < 0
        n0 = length(pts)
        reverse!(pts)
        edge_map = [edge_map[mod1(n0 - i, n0)] for i in 1:n0]
    end
    fronts = [(copy(pts), copy(edge_map))]
    n_events = 0
    while !isempty(fronts)
        fronts = [(p, e) for (fp, em) in fronts
                  for (p, e) in (_drop_reversals(fp, em),)
                  if length(p) >= 3]
        isempty(fronts) && break
        n_events += 1
        n_events <= max_events || error(
            "straight skeleton: no convergence after $max_events wavefront events")
        motion = map(fronts) do (fp, _)
            _compute_bisectors(fp, length(fp))
        end
        t_min = Inf
        for (f, (fp, _)) in enumerate(fronts)
            bis, spd = motion[f]
            tc, _ = _find_all_collapses(fp, length(fp), bis, spd)
            ts, _ = _split_events(fp, bis, spd)
            t_min = min(t_min, tc, ts)
        end
        if t_min == Inf
            shown = length(fronts[1][1]) <= 12 ?
                [(round(p[1]; digits = 3), round(p[2]; digits = 3)) for p in fronts[1][1]] :
                length(fronts[1][1])
            error("straight skeleton: a wavefront has no further event ($shown)")
        end
        next_fronts = Tuple{Vector{NTuple{2,Float64}}, Vector{Int}}[]
        for (f, (fp, em)) in enumerate(fronts)
            n = length(fp)
            bis, spd = motion[f]
            moved = [(fp[i][1] + bis[i][1] * spd[i] * t_min,
                      fp[i][2] + bis[i][2] * spd[i] * t_min) for i in 1:n]
            for i in 1:n
                _dist(fp[i], moved[i]) < 1e-10 && continue
                _record_skeleton_arc!(dcel, registry, fp[i], moved[i], em[mod1(i - 1, n)], em[i])
            end

            ts, splits = _split_events(fp, bis, spd)
            tol = max(1e-9, 1e-6 * t_min)
            pieces = if !isempty(splits) && ts <= t_min + tol
                i, j, _, hit = splits[1]
                # Vertex i becomes `hit` in both pieces; edge j is in both.
                ring(a, b) = a <= b ? collect(a:b) : vcat(collect(a:n), collect(1:b))
                one = ring(mod1(i + 1, n), j)
                two = ring(mod1(j + 1, n), mod1(i - 1, n))
                [(vcat([hit], moved[one]), vcat([em[i]], em[one])),
                 (vcat([hit], moved[two]), vcat([em[j]], em[two]))]
            else
                tc, collapses = _find_all_collapses(fp, n, bis, spd)
                if isempty(collapses) || tc > t_min + tol
                    [(moved, em)]
                else
                    cset = Set(c[1] for c in collapses)
                    cpt = Dict(c[1] => c[2] for c in collapses)
                    if length(cset) == n
                        Tuple{Vector{NTuple{2,Float64}}, Vector{Int}}[]
                    else
                        np, nem, _ = _build_collapsed_polygon_weighted(moved, em, ones(n), n, cset, cpt)
                        [(np, nem)]
                    end
                end
            end
            for (pp, pe) in pieces
                pp, pe = _drop_tiny_edges(pp, pe)
                if length(pp) >= 3
                    push!(next_fronts, (pp, pe))
                else
                    _record_ridge!(dcel, registry, pp, pe)
                end
            end
        end
        fronts = next_fronts
    end
    return dcel
end


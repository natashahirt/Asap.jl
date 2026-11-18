# A supported hole (a beam ring) is a second wavefront. Its edges move into
# the slab at the same unit speed as the outer boundary, and the straight
# skeleton is where those fronts meet. A free opening is not a wavefront:
# the caller still cuts it out of the area, and it carries nothing.

"""Clockwise hole front. Edge `j` of `loop`, as given, keeps face id `id0 + j`."""
function _cw_hole_front(loop::Vector{NTuple{2,Float64}}, id0::Int)
    nh = length(loop)
    ends = NTuple{2, NTuple{2,Float64}}[(loop[j], loop[mod1(j + 1, nh)]) for j in 1:nh]
    if _polygon_area(loop) > 0
        pts = reverse(loop)
        # reverse(loop) walks each given edge backwards. Edge i of `pts`
        # is given edge mod1(nh - i, nh).
        em = Int[id0 + mod1(nh - i, nh) for i in 1:nh]
    else
        pts = copy(loop)
        em = collect((id0 + 1):(id0 + nh))
    end
    return pts, em, ends
end

"""Boundary halfedges of one clockwise hole. The face sits on the left."""
function _add_hole_boundary!(dcel::DCEL, registry::VertexRegistry, pts, em)
    vs = Int[get_or_create_vertex!(registry, p) for p in pts]
    n = length(pts)
    for i in 1:n
        create_halfedge_pair!(dcel, vs[i], vs[mod1(i + 1, n)], em[i], 0)
    end
    return
end

_rot_front(pts, em, k) = k == 1 ? (pts, em) :
    (vcat(pts[k:end], pts[1:(k - 1)]), vcat(em[k:end], em[1:(k - 1)]))

function _edge_frame(pts, i)
    a = pts[i]
    b = pts[mod1(i + 1, length(pts))]
    dx, dy = b[1] - a[1], b[2] - a[2]
    L = hypot(dx, dy)
    L < 1e-12 && return nothing
    return (a = a, b = b, u = (dx / L, dy / L), n = (-dy / L, dx / L), L = L)
end

"""Time until two anti-parallel edges meet, or `Inf`."""
function _edge_meet_time(ptsA, i, ptsB, j)
    ea = _edge_frame(ptsA, i)
    eb = _edge_frame(ptsB, j)
    (ea === nothing || eb === nothing) && return Inf
    dotn = ea.n[1] * eb.n[1] + ea.n[2] * eb.n[2]
    dotn < -1 + 1e-4 || return Inf
    g1 = (eb.a[1] - ea.a[1]) * ea.n[1] + (eb.a[2] - ea.a[2]) * ea.n[2]
    g2 = (eb.b[1] - ea.a[1]) * ea.n[1] + (eb.b[2] - ea.a[2]) * ea.n[2]
    (g1 > 1e-7 && g2 > 1e-7) || return Inf
    abs(g1 - g2) > 1e-4 * max(g1, 1.0) && return Inf
    return (g1 + g2) / 4
end

function _earliest_edge_meet(fronts)
    t = Inf
    for a in eachindex(fronts), b in (a + 1):length(fronts)
        pa, pb = fronts[a][1], fronts[b][1]
        for i in eachindex(pa), j in eachindex(pb)
            t = min(t, _edge_meet_time(pa, i, pb, j))
        end
    end
    return t
end

"""Overlap of two anti-parallel edges that already occupy the same line."""
function _line_overlap(ptsA, i, ptsB, j; tol::Float64 = 1e-6)
    ea = _edge_frame(ptsA, i)
    eb = _edge_frame(ptsB, j)
    (ea === nothing || eb === nothing) && return nothing
    dotn = ea.n[1] * eb.n[1] + ea.n[2] * eb.n[2]
    dotn < -1 + 1e-4 || return nothing
    g1 = (eb.a[1] - ea.a[1]) * ea.n[1] + (eb.a[2] - ea.a[2]) * ea.n[2]
    g2 = (eb.b[1] - ea.a[1]) * ea.n[1] + (eb.b[2] - ea.a[2]) * ea.n[2]
    max(abs(g1), abs(g2)) > tol && return nothing
    s(p) = (p[1] - ea.a[1]) * ea.u[1] + (p[2] - ea.a[2]) * ea.u[2]
    lo = max(0.0, min(s(eb.a), s(eb.b)))
    hi = min(ea.L, max(s(eb.a), s(eb.b)))
    hi - lo < tol && return nothing
    plo = (ea.a[1] + ea.u[1] * lo, ea.a[2] + ea.u[2] * lo)
    phi = (ea.a[1] + ea.u[1] * hi, ea.a[2] + ea.u[2] * hi)
    contained = lo <= tol && hi >= ea.L - tol
    contains = s(eb.a) >= -tol && s(eb.a) <= ea.L + tol &&
               s(eb.b) >= -tol && s(eb.b) <= ea.L + tol
    return (plo = plo, phi = phi, a_contains_b = contains, b_contains_a = contained)
end

function _coincident_pairs(fronts; tol::Float64 = 1e-6)
    pairs = Tuple{Int,Int,Int,Int}[]
    used = Set{Tuple{Int,Int}}()
    for a in eachindex(fronts), b in (a + 1):length(fronts)
        pa, pb = fronts[a][1], fronts[b][1]
        for i in eachindex(pa), j in eachindex(pb)
            ((a, i) in used || (b, j) in used) && continue
            _line_overlap(pa, i, pb, j; tol = tol) === nothing && continue
            push!(pairs, (a, i, b, j))
            push!(used, (a, i), (b, j))
        end
    end
    return pairs
end

function _covers_all(fronts, pairs)
    used = [Set{Int}() for _ in fronts]
    for (a, i, b, j) in pairs
        push!(used[a], i)
        push!(used[b], j)
    end
    return all(f -> length(used[f]) == length(fronts[f][1]), eachindex(fronts))
end

"""
Record the contact and drop both fronts when every edge has met its partner.
Walking against edge A's direction leaves A's face on the left.
"""
function _consume_pairs!(dcel, registry, fronts, pairs)
    for (a, i, b, j) in pairs
        pa, ea = fronts[a]
        a1, a2 = pa[i], pa[mod1(i + 1, length(pa))]
        _record_skeleton_arc!(dcel, registry, a2, a1, ea[i], fronts[b][2][j])
    end
    return Tuple{Vector{NTuple{2,Float64}}, Vector{Int}}[]
end

"""
Splice the shorter edge's front into the longer edge's front along their
overlap. The shorter edge is removed; the longer edge keeps a stub on each
side that still has length, with the same face id.
"""
function _stitch_overlap(longer, il, shorter, is_, dcel, registry)
    pa, ea = _rot_front(longer[1], longer[2], il)
    pb, eb = _rot_front(shorter[1], shorter[2], is_)
    ov = _line_overlap(pa, 1, pb, 1)
    ov === nothing && error("beam-ring wavefronts were paired where their edges do not meet")
    ov.a_contains_b || error(
        "a beam ring met a slab edge in a partial overlap the skeleton does not stitch")
    # Against A's direction: A's face on the left, B's face on the right.
    _record_skeleton_arc!(dcel, registry, ov.phi, ov.plo, ea[1], eb[1])

    s(p) = begin
        ea1 = _edge_frame(pa, 1)
        (p[1] - ea1.a[1]) * ea1.u[1] + (p[2] - ea1.a[2]) * ea1.u[2]
    end
    at_lo_is_b2 = abs(s(pb[2]) - s(ov.plo)) <= abs(s(pb[1]) - s(ov.plo))
    chain = vcat(pb[2:end], pb[1:1])
    chain_em = eb[2:end]
    if !at_lo_is_b2
        reverse!(chain)
        reverse!(chain_em)
    end

    segs = Tuple{NTuple{2,Float64}, NTuple{2,Float64}, Int}[]
    add(p, q, f) = _dist(p, q) >= 1e-8 && push!(segs, (p, q, f))
    add(pa[1], ov.plo, ea[1])
    for k in 1:(length(chain) - 1)
        add(chain[k], chain[k + 1], chain_em[k])
    end
    add(ov.phi, pa[2], ea[1])
    for k in 2:(length(pa) - 1)
        add(pa[k], pa[k + 1], ea[k])
    end
    add(pa[end], pa[1], ea[end])
    isempty(segs) && return (NTuple{2,Float64}[], Int[])

    pts = NTuple{2,Float64}[segs[1][1]]
    em = Int[]
    for (p, q, f) in segs
        _dist(pts[end], p) > 1e-5 && error(
            "beam-ring stitch left a gap in the wavefront")
        push!(em, f)
        _dist(q, pts[1]) > 1e-8 && push!(pts, q)
    end
    length(pts) == length(em) || error(
        "beam-ring stitch produced $(length(pts)) vertices and $(length(em)) edges")
    return _drop_tiny_edges(pts, em)
end

function _move_fronts!(dcel, registry, fronts, motions, t)
    out = Tuple{Vector{NTuple{2,Float64}}, Vector{Int}}[]
    for (f, (pts, em)) in enumerate(fronts)
        bis, spd = motions[f]
        n = length(pts)
        moved = NTuple{2,Float64}[(pts[i][1] + bis[i][1] * spd[i] * t,
                                   pts[i][2] + bis[i][2] * spd[i] * t) for i in 1:n]
        for i in 1:n
            _dist(pts[i], moved[i]) < 1e-10 && continue
            _record_skeleton_arc!(dcel, registry, pts[i], moved[i],
                                  em[mod1(i - 1, n)], em[i])
        end
        push!(out, (moved, em))
    end
    return out
end

"""A vertex of one front hitting the interior of another front's edge."""
function _earliest_vertex_hit(fronts, motions)
    best_t = Inf
    best = (0, 0, 0, 0, (0.0, 0.0))
    for a in eachindex(fronts)
        pa = fronts[a][1]
        bis, spd = motions[a]
        for i in eachindex(pa)
            spd[i] > 0 || continue
            w = bis[i] .* spd[i]
            p = pa[i]
            for b in eachindex(fronts)
                a == b && continue
                pb = fronts[b][1]
                bisb, spdb = motions[b]
                for j in eachindex(pb)
                    k = mod1(j + 1, length(pb))
                    ea, ebpt = pb[j], pb[k]
                    L = _dist(ea, ebpt)
                    L < 1e-12 && continue
                    u = ((ebpt[1] - ea[1]) / L, (ebpt[2] - ea[2]) / L)
                    nrm = (-u[2], u[1])
                    d0 = (p[1] - ea[1]) * nrm[1] + (p[2] - ea[2]) * nrm[2]
                    den = 1.0 - (w[1] * nrm[1] + w[2] * nrm[2])
                    (d0 > 1e-8 && den > 1e-8) || continue
                    t = d0 / den
                    (t > 1e-10 && t < best_t) || continue
                    hit = (p[1] + w[1] * t, p[2] + w[2] * t)
                    at = pb[j] .+ (bisb[j] .* (spdb[j] * t))
                    bt = pb[k] .+ (bisb[k] .* (spdb[k] * t))
                    span = (bt[1] - at[1]) * u[1] + (bt[2] - at[2]) * u[2]
                    s = (hit[1] - at[1]) * u[1] + (hit[2] - at[2]) * u[2]
                    tol = 1e-7 * max(L, 1.0)
                    (span > tol && tol < s < span - tol) || continue
                    best_t = t
                    best = (a, i, b, j, hit)
                end
            end
        end
    end
    return (t = best_t, a = best[1], ia = best[2], b = best[3], jb = best[4], hit = best[5])
end

"""Insert front A's vertex into the interior of front B's edge, one wavefront."""
function _merge_vertex_hit(A, ia, B, jb, hit)
    pb, eb = _rot_front(B[1], B[2], jb)
    pa, ea = _rot_front(A[1], A[2], ia)
    na, nb = length(pa), length(pb)
    segs = Tuple{NTuple{2,Float64}, NTuple{2,Float64}, Int}[]
    add(p, q, f) = _dist(p, q) >= 1e-8 && push!(segs, (p, q, f))
    add(pb[1], hit, eb[1])
    add(hit, pa[2], ea[1])
    for k in 2:(na - 1)
        add(pa[k], pa[k + 1], ea[k])
    end
    add(pa[na], hit, ea[na])
    add(hit, pb[2], eb[1])
    for k in 2:(nb - 1)
        add(pb[k], pb[k + 1], eb[k])
    end
    add(pb[nb], pb[1], eb[nb])
    pts = NTuple{2,Float64}[segs[1][1]]
    em = Int[]
    for (p, q, f) in segs
        _dist(pts[end], p) > 1e-5 && error("beam-ring vertex hit left a gap in the wavefront")
        push!(em, f)
        _dist(q, pts[1]) > 1e-8 && push!(pts, q)
    end
    return _drop_tiny_edges(pts, em)
end

function _earliest_intra(fronts, motions)
    t = Inf
    for (f, (pts, _)) in enumerate(fronts)
        length(pts) < 3 && continue
        bis, spd = motions[f]
        tc, _ = _find_all_collapses(pts, length(pts), bis, spd)
        ts, _ = _split_events(pts, bis, spd)
        t = min(t, tc, ts)
    end
    return t
end

"""Move every front to `t` and apply the intra-front event on each one."""
function _apply_intra(dcel, registry, fronts, motions, t)
    out = Tuple{Vector{NTuple{2,Float64}}, Vector{Int}}[]
    for (f, (fp, em)) in enumerate(fronts)
        n = length(fp)
        bis, spd = motions[f]
        moved = NTuple{2,Float64}[(fp[i][1] + bis[i][1] * spd[i] * t,
                                   fp[i][2] + bis[i][2] * spd[i] * t) for i in 1:n]
        for i in 1:n
            _dist(fp[i], moved[i]) < 1e-10 && continue
            _record_skeleton_arc!(dcel, registry, fp[i], moved[i], em[mod1(i - 1, n)], em[i])
        end
        ts, splits = _split_events(fp, bis, spd)
        tol = max(1e-9, 1e-6 * t)
        pieces = if !isempty(splits) && ts <= t + tol
            i, j, _, hit = splits[1]
            ring(a, b) = a <= b ? collect(a:b) : vcat(collect(a:n), collect(1:b))
            one = ring(mod1(i + 1, n), j)
            two = ring(mod1(j + 1, n), mod1(i - 1, n))
            [(vcat([hit], moved[one]), vcat([em[i]], em[one])),
             (vcat([hit], moved[two]), vcat([em[j]], em[two]))]
        else
            tc, collapses = _find_all_collapses(fp, n, bis, spd)
            if isempty(collapses) || tc > t + tol
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
                push!(out, (pp, pe))
            else
                _record_ridge!(dcel, registry, pp, pe)
            end
        end
    end
    return out
end

"""
Cut front edge `i` at `points` that lie on it.

Each piece keeps the edge's face id. A point on an endpoint is not a cut.
Returns the same front when nothing was inserted.
"""
function _split_front_edge(front, i::Int, points)
    pts, em = front
    n = length(pts)
    ea = _edge_frame(pts, i)
    ea === nothing && return front
    fracs = Float64[]
    for p in points
        t = ((p[1] - ea.a[1]) * ea.u[1] + (p[2] - ea.a[2]) * ea.u[2]) / ea.L
        1e-6 < t < 1 - 1e-6 && push!(fracs, t)
    end
    sort!(fracs)
    unique!(t -> round(t; digits = 8), fracs)
    isempty(fracs) && return front
    a, b = pts[i], pts[mod1(i + 1, n)]
    mids = [(a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])) for t in fracs]
    new_pts = NTuple{2, Float64}[]
    new_em = Int[]
    for k in 1:n
        push!(new_pts, pts[k])
        if k == i
            append!(new_pts, mids)
            for _ in 1:(length(mids) + 1)
                push!(new_em, em[i])
            end
        else
            push!(new_em, em[k])
        end
    end
    return (new_pts, new_em)
end

"""
Shrink the outer front and each hole front together.

An edge of one front meeting an edge of another is removed and the fronts
are joined along that contact. A vertex hitting another front's edge joins
the fronts at that point. With one front left, the ordinary skeleton
finishes it.
"""
function _propagate_with_holes!(dcel::DCEL, registry::VertexRegistry,
                                outer_pts::Vector{NTuple{2,Float64}}, outer_em::Vector{Int},
                                hole_fronts; max_events::Int = 200)
    fronts = Tuple{Vector{NTuple{2,Float64}}, Vector{Int}}[(copy(outer_pts), copy(outer_em))]
    for (pts, em) in hole_fronts
        push!(fronts, (copy(pts), copy(em)))
    end
    n_events = 0
    while !isempty(fronts)
        n_events += 1
        n_events <= max_events || error(
            "straight skeleton: a beam ring did not meet the slab after $max_events events")
        motions = [_compute_bisectors(fp, length(fp)) for (fp, _) in fronts]
        t_edge = _earliest_edge_meet(fronts)
        hit = _earliest_vertex_hit(fronts, motions)
        t_in = _earliest_intra(fronts, motions)
        t = min(t_edge, hit.t, t_in)
        isfinite(t) || error("straight skeleton: a beam ring's wavefront has no further event")
        if t_edge <= hit.t + 1e-8 && t_edge <= t_in + 1e-8
            fronts = _move_fronts!(dcel, registry, fronts, motions, t_edge)
            pairs = _coincident_pairs(fronts)
            isempty(pairs) && error(
                "straight skeleton: opposing edges were due to meet and did not")
            if _covers_all(fronts, pairs)
                fronts = _consume_pairs!(dcel, registry, fronts, pairs)
            else
                a, i, b, j = pairs[1]
                ov_ab = _line_overlap(fronts[a][1], i, fronts[b][1], j)
                # Each edge sticks out past the other. Cut both at the overlap
                # so the contained stitch meets two equal edges. The fronts
                # are already on the contact line, so the pair is recomputed
                # here rather than on the next event.
                if ov_ab !== nothing && !ov_ab.a_contains_b && !ov_ab.b_contains_a
                    fa = _split_front_edge(fronts[a], i, (ov_ab.plo, ov_ab.phi))
                    fb = _split_front_edge(fronts[b], j, (ov_ab.plo, ov_ab.phi))
                    if fa === fronts[a] && fb === fronts[b]
                        error("a beam ring met a slab edge in a partial overlap the skeleton does not stitch")
                    end
                    fronts[a] = fa
                    fronts[b] = fb
                    pairs = _coincident_pairs(fronts)
                    isempty(pairs) && error(
                        "straight skeleton: splitting a partial overlap left no coincident edges")
                    a, i, b, j = pairs[1]
                    ov_ab = _line_overlap(fronts[a][1], i, fronts[b][1], j)
                end
                longer, il, shorter, is_ = ov_ab.a_contains_b ?
                    (fronts[a], i, fronts[b], j) : (fronts[b], j, fronts[a], i)
                stitched = _stitch_overlap(longer, il, shorter, is_, dcel, registry)
                keep = [fronts[k] for k in eachindex(fronts) if k != a && k != b]
                length(stitched[1]) >= 3 && push!(keep, stitched)
                fronts = keep
            end
        elseif hit.t <= t_in + 1e-8
            fronts = _move_fronts!(dcel, registry, fronts, motions, hit.t)
            # The hit was computed on the unmoved front. Re-find it on the
            # moved front: the vertex now sits on the edge.
            hit2 = _vertex_on_edge(fronts)
            hit2 === nothing && error("straight skeleton: a ring vertex missed the slab edge")
            a, ia, b, jb, pt = hit2
            merged = _merge_vertex_hit(fronts[a], ia, fronts[b], jb, pt)
            keep = [fronts[k] for k in eachindex(fronts) if k != a && k != b]
            length(merged[1]) >= 3 && push!(keep, merged)
            fronts = keep
        else
            fronts = _apply_intra(dcel, registry, fronts, motions, t_in)
        end
        if length(fronts) <= 1
            if length(fronts) == 1 && length(fronts[1][1]) >= 3
                _propagate_isotropic!(dcel, registry, fronts[1][1], fronts[1][2])
            elseif length(fronts) == 1
                _record_ridge!(dcel, registry, fronts[1][1], fronts[1][2])
            end
            return dcel
        end
    end
    return dcel
end

"""A vertex of one front lying on an edge of another, after the move."""
function _vertex_on_edge(fronts)
    for a in eachindex(fronts), b in eachindex(fronts)
        a == b && continue
        pa, pb = fronts[a][1], fronts[b][1]
        for i in eachindex(pa), j in eachindex(pb)
            k = mod1(j + 1, length(pb))
            ea, eb = pb[j], pb[k]
            L = _dist(ea, eb)
            L < 1e-9 && continue
            u = ((eb[1] - ea[1]) / L, (eb[2] - ea[2]) / L)
            nrm = (-u[2], u[1])
            p = pa[i]
            d = (p[1] - ea[1]) * nrm[1] + (p[2] - ea[2]) * nrm[2]
            abs(d) > 1e-6 && continue
            s = (p[1] - ea[1]) * u[1] + (p[2] - ea[2]) * u[2]
            tol = 1e-6 * max(L, 1.0)
            tol < s < L - tol || continue
            return (a, i, b, j, p)
        end
    end
    return nothing
end

"""Tributary polygons of the ring edges, in the order the loops were given."""
function _ring_tributaries(dcel::DCEL, pts_orig, faces::Vector{Int},
                           ends::Vector{NTuple{2, NTuple{2,Float64}}})
    out = TributaryPolygon[]
    m = length(pts_orig)
    for k in eachindex(faces)
        poly = extract_face_polygon(dcel, faces[k])
        length(poly) < 3 && (poly = extract_face_polygon_by_edges(dcel, faces[k]))
        area = length(poly) < 3 ? 0.0 : abs(_polygon_area(poly))
        push!(out, _make_tributary_with_depths(m + k, poly, pts_orig, area, 0.0; ends = ends[k]))
    end
    return out
end

function _outer_and_ring_lengths(pts_orig, ends)
    n = length(pts_orig)
    L = Float64[hypot(pts_orig[mod1(i + 1, n)][1] - pts_orig[i][1],
                      pts_orig[mod1(i + 1, n)][2] - pts_orig[i][2]) for i in 1:n]
    for (a, b) in ends
        push!(L, hypot(b[1] - a[1], b[2] - a[2]))
    end
    return L
end

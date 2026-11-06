# =============================================================================
# Helper Functions (DCEL-specific versions)
# =============================================================================

"""Distance between two points."""
_dist(p1, p2) = hypot(p1[1] - p2[1], p1[2] - p2[2])

"""
Clip polygon poly by half-plane: dot(x - p0, n) >= 0
Sutherland–Hodgman clipping against a half-plane.
"""
function clip_halfplane(poly::Vector{NTuple{2,Float64}},
                        p0::NTuple{2,Float64},
                        n::NTuple{2,Float64};
                        eps=1e-12)
    isempty(poly) && return poly

    function inside(p)
        ((p[1]-p0[1])*n[1] + (p[2]-p0[2])*n[2]) >= -eps
    end

    function intersect(a, b)
        # find t in [0,1] where dot((a + t(b-a)) - p0, n) == 0
        ax = (a[1]-p0[1])*n[1] + (a[2]-p0[2])*n[2]
        bx = (b[1]-p0[1])*n[1] + (b[2]-p0[2])*n[2]
        denom = bx - ax
        if abs(denom) < eps
            return b  # nearly parallel; fallback
        end
        t = -ax / denom
        t = clamp(t, 0.0, 1.0)
        return (a[1] + t*(b[1]-a[1]), a[2] + t*(b[2]-a[2]))
    end

    out = NTuple{2,Float64}[]
    nV = length(poly)
    for i in 1:nV
        a = poly[i]
        b = poly[mod1(i+1, nV)]
        ina = inside(a)
        inb = inside(b)
        if ina && inb
            push!(out, b)
        elseif ina && !inb
            push!(out, intersect(a,b))
        elseif !ina && inb
            push!(out, intersect(a,b))
            push!(out, b)
        end
    end

    return out
end

"""
Split a merged-face polygon P for a collinear chain p0..pk into per-subedge polygons.
Returns a vector of polygons, one for each sub-edge in the chain.
"""
function split_collinear_face(P::Vector{NTuple{2,Float64}},
                              chain_pts::Vector{NTuple{2,Float64}})
    k = length(chain_pts) - 1
    k >= 1 || return [P]

    # tangent along the chain
    tvec = (chain_pts[end][1]-chain_pts[1][1], chain_pts[end][2]-chain_pts[1][2])
    tl = hypot(tvec...)
    tl < 1e-15 && error("degenerate chain")
    t = (tvec[1]/tl, tvec[2]/tl)

    polys = Vector{Vector{NTuple{2,Float64}}}(undef, k)
    for i in 1:k
        Pi = P
        # keep points with dot(x - p_i, t) >= 0
        Pi = clip_halfplane(Pi, chain_pts[i], t)
        # keep points with dot(x - p_{i+1}, t) <= 0  <=> dot(x - p_{i+1}, -t) >= 0
        Pi = clip_halfplane(Pi, chain_pts[i+1], (-t[1], -t[2]))
        polys[i] = Pi
    end
    return polys
end


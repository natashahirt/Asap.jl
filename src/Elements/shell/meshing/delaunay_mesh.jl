#=
Shell Creation Utilities
========================

Convenient functions to create shell elements from corner nodes or external meshes.
Uses DelaunayTriangulation.jl for automatic polygon triangulation.
Includes structured mesh generation (T3block) adapted from FinEtools.jl (Petr Krysl, MIT License).
=#

import DelaunayTriangulation as DT
import Meshes
import Meshes: SimpleMesh, Point, coords
import Random: MersenneTwister

# DelaunayTriangulation draws from an RNG during constrained insertion and
# Ruppert refinement. A fresh generator with this seed is passed on every
# call, so a mesh is a function of its points and does not depend on
# `Random.default_rng()` or on meshes built earlier in the process.
const _ASAP_MESH_SEED = 0x41534150
const _MESH_RNG_CALLS = Ref(0)
function _mesh_rng()
    _MESH_RNG_CALLS[] += 1
    return MersenneTwister(_ASAP_MESH_SEED)
end

"""
Exact nearest-neighbor distance. Matches

    minimum(hypot(x - p[1], y - p[2]) for p in pts)

The grid only decides which points to measure. Search starts at the first
ring that can hold a point, then grows until every unsearched cell is at
least as far as the best distance already found. A tight cluster keeps a
cell of 1 so a far query does not walk empty rings.
"""
struct _PointGrid
    originx::Float64
    originy::Float64
    cell::Float64
    buckets::Dict{NTuple{2, Int}, Vector{Int}}
    pts::Vector{Tuple{Float64, Float64}}
    imin::Int
    imax::Int
    jmin::Int
    jmax::Int
end

function _PointGrid(pts::Vector{Tuple{Float64, Float64}})
    n = length(pts)
    empty_buckets = Dict{NTuple{2, Int}, Vector{Int}}()
    n == 0 && return _PointGrid(0.0, 0.0, 1.0, empty_buckets, pts, 0, 0, 0, 0)
    xmin = ymin = Inf
    xmax = ymax = -Inf
    for (x, y) in pts
        xmin = min(xmin, x); xmax = max(xmax, x)
        ymin = min(ymin, y); ymax = max(ymax, y)
    end
    span = max(xmax - xmin, ymax - ymin)
    # A coincident cluster (one column, or targets closer than 0.1 mm) uses a
    # 1 m cell. A cell built from a near-zero span makes the ring index enormous
    # and the search never reaches the point.
    cell = span > 1e-4 ? span / max(sqrt(float(n)), 1.0) : 1.0
    buckets = empty_buckets
    imin = imax = jmin = jmax = 0
    for (i, (x, y)) in enumerate(pts)
        key = (floor(Int, (x - xmin) / cell), floor(Int, (y - ymin) / cell))
        push!(get!(() -> Int[], buckets, key), i)
        if i == 1
            imin = imax = key[1]
            jmin = jmax = key[2]
        else
            imin = min(imin, key[1]); imax = max(imax, key[1])
            jmin = min(jmin, key[2]); jmax = max(jmax, key[2])
        end
    end
    return _PointGrid(xmin, ymin, cell, buckets, pts, imin, imax, jmin, jmax)
end

function _consider_bucket!(g::_PointGrid, ix::Int, iy::Int, x::Float64, y::Float64, best::Float64)
    bucket = get(g.buckets, (ix, iy), nothing)
    bucket === nothing && return best
    for i in bucket
        px, py = g.pts[i]
        d = hypot(x - px, y - py)
        d < best && (best = d)
    end
    return best
end

function _min_dist(g::_PointGrid, x::Float64, y::Float64)
    isempty(g.pts) && return Inf
    ix = floor(Int, (x - g.originx) / g.cell)
    iy = floor(Int, (y - g.originy) / g.cell)
    dx0 = ix < g.imin ? g.imin - ix : ix > g.imax ? ix - g.imax : 0
    dy0 = iy < g.jmin ? g.jmin - iy : iy > g.jmax ? iy - g.jmax : 0
    ring = max(dx0, dy0)
    far = max(
        max(abs(ix - g.imin), abs(iy - g.jmin)),
        max(abs(ix - g.imin), abs(iy - g.jmax)),
        max(abs(ix - g.imax), abs(iy - g.jmin)),
        max(abs(ix - g.imax), abs(iy - g.jmax)),
    )
    best = Inf
    while ring <= far
        if ring == 0
            best = _consider_bucket!(g, ix, iy, x, y, best)
        else
            y0 = iy - ring
            y1 = iy + ring
            for dx in -ring:ring
                best = _consider_bucket!(g, ix + dx, y0, x, y, best)
                best = _consider_bucket!(g, ix + dx, y1, x, y, best)
            end
            for dy in (-ring + 1):(ring - 1)
                best = _consider_bucket!(g, ix - ring, iy + dy, x, y, best)
                best = _consider_bucket!(g, ix + ring, iy + dy, x, y, best)
            end
        end
        best <= ring * g.cell && return best
        ring += 1
    end
    return best
end


struct ForceEnvelopes
    element::Element
    resolution::Integer
    x::Vector{Float64}
    Plow::Vector{Float64}
    Phigh::Vector{Float64}
    Mylow::Vector{Float64}
    Myhigh::Vector{Float64}
    Vylow::Vector{Float64}
    Vyhigh::Vector{Float64}
    Mzlow::Vector{Float64}
    Mzhigh::Vector{Float64}
    Vzlow::Vector{Float64}
    Vzhigh::Vector{Float64}
end

"""
    load_envelopes(model, loads::Vector{Vector{<:AbstractLoad}})

Get the high/low internal forces for a series of external loads
"""
function load_envelopes(model::Union{FrameModel, Model}, loads::Vector{Vector{<:AbstractLoad}}, increment::Real)
    envelopes = Vector{ForceEnvelopes}()
    forceresults = Vector{Vector{InternalForces}}()

    for load in loads
        solve!(model, load)
        push!(forceresults, forces(model, increment))
    end

    n_elem = length(first(forceresults))
    n_loads = length(forceresults)
    sizehint!(envelopes, n_elem)

    for i = 1:n_elem
        ref = first(forceresults)[i]
        e   = ref.element
        res = ref.resolution
        x   = ref.x
        n_pts = length(x)

        # Preallocate matrices (pts × load cases) — avoids hcat+getproperty broadcasts
        Pm  = Matrix{Float64}(undef, n_pts, n_loads)
        Mym = Matrix{Float64}(undef, n_pts, n_loads)
        Vym = Matrix{Float64}(undef, n_pts, n_loads)
        Mzm = Matrix{Float64}(undef, n_pts, n_loads)
        Vzm = Matrix{Float64}(undef, n_pts, n_loads)

        @inbounds for j in 1:n_loads
            fr = forceresults[j][i]
            for k in 1:n_pts
                Pm[k, j]  = fr.P[k]
                Mym[k, j] = fr.My[k]
                Vym[k, j] = fr.Vy[k]
                Mzm[k, j] = fr.Mz[k]
                Vzm[k, j] = fr.Vz[k]
            end
        end

        # Compute envelopes row-by-row
        Plow  = Vector{Float64}(undef, n_pts);  Phigh  = Vector{Float64}(undef, n_pts)
        Mylow = Vector{Float64}(undef, n_pts);  Myhigh = Vector{Float64}(undef, n_pts)
        Vylow = Vector{Float64}(undef, n_pts);  Vyhigh = Vector{Float64}(undef, n_pts)
        Mzlow = Vector{Float64}(undef, n_pts);  Mzhigh = Vector{Float64}(undef, n_pts)
        Vzlow = Vector{Float64}(undef, n_pts);  Vzhigh = Vector{Float64}(undef, n_pts)

        @inbounds for k in 1:n_pts
            lo, hi = Pm[k, 1], Pm[k, 1]
            for j in 2:n_loads; v = Pm[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Plow[k] = lo; Phigh[k] = hi

            lo, hi = Mym[k, 1], Mym[k, 1]
            for j in 2:n_loads; v = Mym[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Mylow[k] = lo; Myhigh[k] = hi

            lo, hi = Vym[k, 1], Vym[k, 1]
            for j in 2:n_loads; v = Vym[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Vylow[k] = lo; Vyhigh[k] = hi

            lo, hi = Mzm[k, 1], Mzm[k, 1]
            for j in 2:n_loads; v = Mzm[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Mzlow[k] = lo; Mzhigh[k] = hi

            lo, hi = Vzm[k, 1], Vzm[k, 1]
            for j in 2:n_loads; v = Vzm[k,j]; v < lo && (lo = v); v > hi && (hi = v); end
            Vzlow[k] = lo; Vzhigh[k] = hi
        end

        push!(envelopes, ForceEnvelopes(e, res, x,
            Plow, Phigh, Mylow, Myhigh, Vylow, Vyhigh, Mzlow, Mzhigh, Vzlow, Vzhigh))
    end

    return envelopes
end

# =============================================================================
# Convenience: moments! / shear! for frame elements
# =============================================================================

"""
    moments!(My, Mz, element, model; resolution=20)

Write bending moments My (about local Y) and Mz (about local Z) into
pre-allocated vectors.  Returns `(My, Mz)`.
"""
function moments!(My::Vector{Float64}, Mz::Vector{Float64},
                  element::AbstractElement, model::Union{FrameModel, Model};
                  resolution::Int = 20)
    eif = ElementInternalForces(element, model; resolution=resolution)
    copyto!(My, eif.My)
    copyto!(Mz, eif.Mz)
    return (My, Mz)
end

"""
    shears!(Vy, Vz, element, model; resolution=20)

Write shear forces Vy and Vz into pre-allocated vectors.  Returns `(Vy, Vz)`.
"""
function shears!(Vy::Vector{Float64}, Vz::Vector{Float64},
                 element::AbstractElement, model::Union{FrameModel, Model};
                 resolution::Int = 20)
    eif = ElementInternalForces(element, model; resolution=resolution)
    copyto!(Vy, eif.Vy)
    copyto!(Vz, eif.Vz)
    return (Vy, Vz)
end

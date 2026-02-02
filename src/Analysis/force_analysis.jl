
"""
Accumulate the internal forces caused by a TributaryLoad to the current element.
TributaryLoad represents a piecewise linear distributed load based on tributary widths.
Each segment between breakpoints is approximated as a partial uniform load.
"""
function accumulate_force!(load::TributaryLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    element = load.element
    R = element.R[1:3, 1:3]
    L = ustrip(u"m", uconvert(u"m", element.length))
    
    # Get load intensity per unit length at each position: w(s) = pressure * width(s)
    # Direction in global coordinates
    dir = collect(load.direction)
    pressure = ustrip(u"Pa", load.pressure)
    widths_m = [ustrip(u"m", w) for w in load.widths]
    
    # Process each segment between breakpoints
    positions = load.positions
    n_seg = length(positions) - 1
    
    for i in 1:n_seg
        s1, s2 = positions[i], positions[i+1]
        w1, w2 = widths_m[i], widths_m[i+1]
        
        # Skip negligible segments
        (s2 - s1) < 1e-12 && continue
        (w1 + w2) < 1e-12 && continue
        
        # Segment bounds in physical coordinates
        a = s1 * L  # start of loaded segment
        b = s2 * L  # end of loaded segment
        
        # Average intensity for this segment (trapezoidal → uniform approximation)
        w_avg = pressure * (w1 + w2) / 2  # N/m
        
        # Global load vector (typically downward: 0, 0, -1)
        w_global = dir .* w_avg
        
        # Transform to local coordinate system
        w_local = R * w_global
        wx, wy, wz = w_local .* [1, -1, -1]
        
        # Accumulate forces from this partial uniform load segment [a, b]
        for (j, x) in enumerate(xvals)
            # Axial force contribution
            P[j] += _PLine_partial(wx, L, x, a, b)
            
            # Bending in y-z plane
            My[j] += _MLine_partial(element, wy, L, x, a, b)
            Vy[j] += _VLine_partial(element, wy, L, x, a, b)
            
            # Bending in x-z plane  
            Mz[j] += _MLine_partial(element, wz, L, x, a, b)
            Vz[j] += _VLine_partial(element, wz, L, x, a, b)
        end
    end
end

# Helper: axial force from partial uniform load on [a, b]
function _PLine_partial(w, L, x, a, b)
    # Total force from segment
    F = w * (b - a)
    # Reaction at start (assuming simply supported for axial)
    R_start = F * (L - (a + b) / 2) / L
    
    if x < a
        -R_start
    elseif x < b
        -(R_start + w * (x - a))
    else
        -(R_start + F)
    end
end

# Helper: moment from partial uniform load on [a, b]
function _MLine_partial(::Element{FreeFree}, w, L, x, a, b)
    _MLine_partial_ss(w, L, x, a, b)
end
function _MLine_partial(::Element{Joist}, w, L, x, a, b)
    _MLine_partial_ss(w, L, x, a, b)
end
function _MLine_partial(::Element{FixedFixed}, w, L, x, a, b)
    _MLine_partial_ff(w, L, x, a, b)
end
function _MLine_partial(::Element{FreeFixed}, w, L, x, a, b)
    _MLine_partial_pf(w, L, x, a, b)
end
function _MLine_partial(::Element{FixedFree}, w, L, x, a, b)
    _MLine_partial_fp(w, L, x, a, b)
end

# Simply supported partial uniform load moment
function _MLine_partial_ss(w, L, x, a, b)
    # Reaction at A (x=0)
    c = (a + b) / 2  # centroid of load
    F = w * (b - a)
    Ra = F * (L - c) / L
    
    if x <= a
        Ra * x
    elseif x <= b
        Ra * x - w * (x - a)^2 / 2
    else
        Ra * x - F * (x - c)
    end
end

# Fixed-fixed partial uniform load moment
function _MLine_partial_ff(w, L, x, a, b)
    seg_len = b - a
    c = (a + b) / 2
    F = w * seg_len
    
    # Fixed-end moments (approximate)
    Ma = F * c * (L - c)^2 / L^2 - F * seg_len^2 / 12 * (L - 2*c) / L
    Mb = F * (L - c) * c^2 / L^2 + F * seg_len^2 / 12 * (L - 2*c) / L
    
    # Reaction at A
    Ra = (F * (L - c) + Ma - Mb) / L
    
    if x <= a
        Ra * x - Ma
    elseif x <= b
        Ra * x - Ma - w * (x - a)^2 / 2
    else
        Ra * x - Ma - F * (x - c)
    end
end

# Pinned-fixed (approximate)
function _MLine_partial_pf(w, L, x, a, b)
    _MLine_partial_ss(w, L, x, a, b)
end

# Fixed-pinned (approximate)
function _MLine_partial_fp(w, L, x, a, b)
    _MLine_partial_ss(w, L, x, a, b)
end

# Helper: shear from partial uniform load on [a, b]
function _VLine_partial(element, w, L, x, a, b)
    c = (a + b) / 2
    F = w * (b - a)
    Ra = F * (L - c) / L  # Simply supported reaction
    
    if x <= a
        Ra
    elseif x <= b
        Ra - w * (x - a)
    else
        Ra - F
    end
end

"""
Accumlate the internal forces cause by a given load to the current element
"""
function accumulate_force!(load::LineLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", load.element.length))

    # distributed load magnitudes in LCS (strip units from load.value)
    value_stripped = [ustrip(u"N/m", uconvert(u"N/m", v)) for v in load.value]
    wx, wy, wz = (R * value_stripped) .* [1, -1, -1]

    P .+= [PLine(wx, L, x) for x in xvals]

    My .+= [MLine(load.element, wy, L, x) for x in xvals]
    Vy .+= [VLine(load.element, wy, L, x) for x in xvals]

    Mz .+= [MLine(load.element, wz, L, x) for x in xvals]
    Vz .+= [VLine(load.element, wz, L, x) for x in xvals]
end

"""
Accumlate the internal forces cause by a given load to the current element
"""
function accumulate_force!(load::PointLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", load.element.length))
    frac = load.position

    # distributed load magnitudes in LCS (strip units from load.value)
    value_stripped = [ustrip(u"N", uconvert(u"N", v)) for v in load.value]
    px, py, pz = (R * value_stripped) .* [1, -1, -1]

    P .+= [PPoint(px, L, x, frac) for x in xvals]

    My .+= [MPoint(load.element, py, L, x, frac) for x in xvals]
    Vy .+= [VPoint(load.element, py, L, x, frac) for x in xvals]

    Mz .+= [MPoint(load.element, pz, L, x, frac) for x in xvals]
    Vz .+= [VPoint(load.element, pz, L, x, frac) for x in xvals]
end

"""
Accumulate the internal forces caused by a GravityLoad (self-weight) to the current element.
GravityLoad applies a distributed load equal to ρ * A * g in the global -Z direction.
"""
function accumulate_force!(load::GravityLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    element = load.element
    R = element.R[1:3, 1:3]
    # Strip units from length
    L = ustrip(u"m", uconvert(u"m", element.length))

    # Self-weight per unit length: w = ρ * A * g (in N/m)
    ρ = ustrip(u"kg/m^3", element.section.ρ)
    A = ustrip(u"m^2", element.section.A)
    g = ustrip(u"m/s^2", load.factor)
    w_mag = ρ * A * g  # N/m

    # Global load vector: self-weight acts in -Z direction
    w_global = [0.0, 0.0, -w_mag]

    # Transform to local coordinate system
    wx, wy, wz = (R * w_global) .* [1, -1, -1]

    P .+= [PLine(wx, L, x) for x in xvals]

    My .+= [MLine(element, wy, L, x) for x in xvals]
    Vy .+= [VLine(element, wy, L, x) for x in xvals]

    Mz .+= [MLine(element, wz, L, x) for x in xvals]
    Vz .+= [VLine(element, wz, L, x) for x in xvals]
end

"""
    ElementInternalForces

Internal force results for a frame element (beam/column).

# Fields
- `element::Element` - Reference to the frame element
- `resolution::Integer` - Number of sampling points
- `x::Vector{Float64}` - Position along element [m]
- `P::Vector{Float64}` - Axial force [N]
- `My::Vector{Float64}` - Bending moment about local Y [N·m]
- `Vy::Vector{Float64}` - Shear force in local Y [N]
- `Mz::Vector{Float64}` - Bending moment about local Z [N·m]
- `Vz::Vector{Float64}` - Shear force in local Z [N]
"""
struct ElementInternalForces
    element::Element
    resolution::Integer
    x::Vector{Float64}
    P::Vector{Float64}
    My::Vector{Float64}
    Vy::Vector{Float64}
    Mz::Vector{Float64}
    Vz::Vector{Float64}
end

# Note: InternalForces is now a unified dispatch function defined in shell_forces.jl
# that returns ElementInternalForces for frame elements and ShellInternalForces for shells.


"""
    ElementInternalForces(element, model; resolution=20)

Internal force sampling for a frame element from a solved model.
"""
function ElementInternalForces(element::AbstractElement, model::Union{FrameModel, Model}; resolution = 20)
    
    #beam information
    dofs = etype2DOF[typeof(element)]
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", element.length))

    #discretization
    xinc = collect(range(0, L, resolution))

    #end node information
    # Strip units from displacements for matrix multiplication (R and K are Float64)
    uglobal = [to_displacement_vec(element.nodeStart.displacement); to_displacement_vec(element.nodeEnd.displacement)]
    
    # end forces that are relevant to the given element/release condition
    Flocal = (element.R * element.K * uglobal) .* dofs

    # shear/moment acting at the *starting* point of an element in LCS
    Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

    # initialize internal force vectors
    P = repeat([Pstart], resolution)
    My = Vystart .* xinc .- Mystart
    Vy = zero(My) .+ Vystart
    Mz = Vzstart .* xinc .- Mzstart
    Vz = zero(Mz) .+ Vzstart

    element_loads = get_elemental_loads(model)

    # accumulate loads
    for load in element_loads[element.elementID]
        accumulate_force!(load,
            xinc,
            P,
            My,
            Vy,
            Mz,
            Vz)  
    end

    return ElementInternalForces(element, resolution, xinc, P, My, Vy, Mz, Vz)
end

"""
    ElementInternalForces(element::Element, loads::Vector{<:ElementLoad}; resolution = 20)

Get internal force results for a given element from a set of loads
"""
function ElementInternalForces(element::AbstractElement, loads::AbstractVector{<:AbstractLoad}; resolution = 20)
    
    #beam information
    dofs = etype2DOF[typeof(element)]
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", element.length))

    #discretization
    xinc = collect(range(0, L, resolution))

    #end node information
    # Strip units from displacements for matrix multiplication (R and K are Float64)
    uglobal = [to_displacement_vec(element.nodeStart.displacement); to_displacement_vec(element.nodeEnd.displacement)]
    
    # end forces that are relevant to the given element/release condition
    Flocal = (element.R * element.K * uglobal) .* dofs

    # shear/moment acting at the *starting* point of an element in LCS
    Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

    # initialize internal force vectors
    P = repeat([Pstart], resolution)
    My = Vystart .* xinc .- Mystart
    Vy = zero(My) .+ Vystart
    Mz = Vzstart .* xinc .- Mzstart
    Vz = zero(Mz) .+ Vzstart

    # accumulate loads
    for load in loads
        accumulate_force!(load,
            xinc,
            P,
            My,
            Vy,
            Mz,
            Vz)  
    end

    return ElementInternalForces(element, resolution, xinc, P, My, Vy, Mz, Vz)
end


function get_elemental_loads(model::FrameModel)
    element_to_loads = [AbstractLoad[] for _ in 1:model.nElements]
    for load in model.loads
        hasproperty(load, :element) || continue
        push!(element_to_loads[load.element.elementID], load)
    end
    return element_to_loads
end

function get_elemental_loads(model::Model)
    element_to_loads = [AbstractLoad[] for _ in 1:model.nFrameElements]
    for load in model.loads
        hasproperty(load, :element) || continue
        # Only include frame element loads
        if load.element isa FrameElement
            push!(element_to_loads[load.element.elementID], load)
        end
    end
    return element_to_loads
end

"""
    ElementInternalForces(elements::Vector{<:FrameElement}, model; resolution=20)

Get internal force results for a group of ordered elements that form a single physical element.
"""
function ElementInternalForces(elements::AbstractVector{<:AbstractElement}, model::Union{FrameModel, Model}; resolution = 20)
    
    xstore = Vector{Float64}()
    pstore = Vector{Float64}()
    mystore = Vector{Float64}()
    vystore = Vector{Float64}()
    mzstore = Vector{Float64}()
    vzstore = Vector{Float64}()

    resolution = Int(round(resolution / length(elements)))

    element_loads_map = get_elemental_loads(model)

    #beam information
    for element in elements
        loadids = element_loads_map[element.elementID]

        dofs = etype2DOF[typeof(element)]
        # Strip units from length (convert Quantity → Float64)
        L = ustrip(u"m", uconvert(u"m", element.length))

        #discretization
        xinc = collect(range(0, L, max(resolution, 2)))

        #end node information
        uglobal = [element.nodeStart.displacement; element.nodeEnd.displacement]
        
        # end forces that are relevant to the given element/release condition
        Flocal = (element.R * element.K * uglobal) .* dofs

        # shear/moment acting at the *starting* point of an element in LCS
        Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

        # initialize internal force vectors
        P = repeat([Pstart], resolution)
        My = Vystart .* xinc .- Mystart
        Vy = zero(My) .+ Vystart
        Mz = Vzstart .* xinc .- Mzstart
        Vz = zero(Mz) .+ Vzstart

        # accumulate loads
        for load in loadids
            accumulate_force!(load,
                xinc,
                P,
                My,
                Vy,
                Mz,
                Vz)  
        end

        if isempty(xstore)
            xstore = [xstore; xinc]
        else
            xstore = [xstore; xstore[end] .+ xinc]
        end

        pstore = [pstore; P]
        mystore = [mystore; My]
        vystore = [vystore; Vy]
        mzstore = [mzstore; Mz]
        vzstore = [vzstore; Vz]

    end

    return ElementInternalForces(elements[1], resolution, xstore, pstore, mystore, vystore, mzstore, vzstore)
end



"""
    ElementInternalForces(model::FrameModel, increment)

Get the internal forces of all frame elements in a FrameModel.

# Arguments
- `model` - FrameModel (frame elements only)
- `increment` - Distance between sampling points (accepts Unitful length or Real in meters)
"""
function ElementInternalForces(model::FrameModel, increment)
    # Accept Unitful or Real - strip to meters internally
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)
    
    element_loads = get_elemental_loads(model)

    #get the internal forces
    IF = [ElementInternalForces(e, model.u[e.globalID], element_loads[e.elementID], inc_m, etype2DOF[typeof(e)]) for e in model.elements]

    #concatenate the results
    x = vcat([IF[i].x .+ IF[i].L / 2 .+ sum(IF[j].L for j = 1:i-1; init = 0.0) for i = 1:length(IF)]...)
    L = sum([IF[i].L for i = 1:length(IF)])
    P = vcat([IF[i].P for i = 1:length(IF)]...)
    My = vcat([IF[i].My for i = 1:length(IF)]...)
    Vy = vcat([IF[i].Vy for i = 1:length(IF)]...)
    Mz = vcat([IF[i].Mz for i = 1:length(IF)]...)
    Vz = vcat([IF[i].Vz for i = 1:length(IF)]...)

    return ElementInternalForces(nothing, L, x, P, My, Vy, Mz, Vz)
end

"""
    ElementInternalForces(model::Model, increment)

Get the internal forces of all frame elements in a Model.

# Arguments
- `model` - Model (may contain both frame and shell elements)
- `increment` - Distance between sampling points (accepts Unitful length or Real in meters)

# Returns
Vector of `ElementInternalForces`, one per frame element.
"""
function ElementInternalForces(model::Model, increment)
    # Accept Unitful or Real - strip to meters internally
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)

    if isempty(model.frame_elements)
        return ElementInternalForces[]
    end

    results = Vector{ElementInternalForces}()
    element_loads_map = get_elemental_loads(model)

    for element in model.frame_elements
        L = ustrip(u"m", element.length)
        n = max(Int(round(L / inc_m)), 2)

        push!(results, ElementInternalForces(element, element_loads_map[element.elementID]; resolution = n))
    end

    return results
end

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
    #
    envelopes = Vector{ForceEnvelopes}()

    #collector of force results
    forceresults = Vector{Vector{InternalForces}}()

    # perform analysis
    for load in loads
        solve!(model, load)
        push!(forceresults, forces(model, increment))
    end

    # number of actual elements
    n = length(first(forceresults))

    # for each element
    for i = 1:n
        e = first(forceresults)[i].element
        res = first(forceresults)[i].resolution
        x = first(forceresults)[i].x

        P = hcat(getproperty.(getindex.(forceresults, i), :P)...)
        My = hcat(getproperty.(getindex.(forceresults, i), :My)...)
        Vy = hcat(getproperty.(getindex.(forceresults, i), :Vy)...)
        Mz = hcat(getproperty.(getindex.(forceresults, i), :Mz)...)
        Vz = hcat(getproperty.(getindex.(forceresults, i), :Vz)...)

        Prange = extrema.(eachrow(P))
        Myrange = extrema.(eachrow(My))
        Vyrange = extrema.(eachrow(Vy))
        Mzrange = extrema.(eachrow(Mz))
        Vzrange = extrema.(eachrow(Vz))

        envelope = ForceEnvelopes(e,
            res,
            x,
            getindex.(Prange, 1),
            getindex.(Prange, 2),
            getindex.(Myrange, 1),
            getindex.(Myrange, 2),
            getindex.(Vyrange, 1),
            getindex.(Vyrange, 2),
            getindex.(Mzrange, 1),
            getindex.(Mzrange, 2),
            getindex.(Vzrange, 1),
            getindex.(Vzrange, 2))

        push!(envelopes, envelope)
    end

    return envelopes
end
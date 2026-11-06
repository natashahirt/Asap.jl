"""
    ElementInternalForces(elements::Vector{<:FrameElement}, model; resolution=20)

Get internal force results for a group of ordered elements that form a single physical element.
"""
function ElementInternalForces(elements::AbstractVector{<:AbstractElement}, model::Union{FrameModel, Model}; resolution = 20)
    
    per_elem_res = max(Int(round(resolution / length(elements))), 2)
    total_pts = per_elem_res * length(elements)

    xstore  = Vector{Float64}(undef, total_pts)
    pstore  = Vector{Float64}(undef, total_pts)
    mystore = Vector{Float64}(undef, total_pts)
    vystore = Vector{Float64}(undef, total_pts)
    mzstore = Vector{Float64}(undef, total_pts)
    vzstore = Vector{Float64}(undef, total_pts)

    # Scratch vectors reused per element
    P_tmp  = Vector{Float64}(undef, per_elem_res)
    My_tmp = Vector{Float64}(undef, per_elem_res)
    Vy_tmp = Vector{Float64}(undef, per_elem_res)
    Mz_tmp = Vector{Float64}(undef, per_elem_res)
    Vz_tmp = Vector{Float64}(undef, per_elem_res)
    xinc   = Vector{Float64}(undef, per_elem_res)

    element_loads_map = get_elemental_loads(model)
    k = 0

    # Preallocated mul! buffers for R*K*u
    _Ku_buf = Vector{Float64}(undef, 12)
    _Fl_buf = Vector{Float64}(undef, 12)

    for element in elements
        loadids = element_loads_map[element.elementID]

        dofs = etype2DOF[typeof(element)]
        ep = ElementProps(element)
        L = ep.L

        # Fill xinc
        rng = range(0, L, per_elem_res)
        @inbounds for j in 1:per_elem_res
            xinc[j] = rng[j]
        end

        # Use model.u directly (Float64) + mul! to avoid intermediates
        u_e = model.u[element.globalID]
        mul!(_Ku_buf, element.K, u_e)
        mul!(_Fl_buf, element.R, _Ku_buf)
        @inbounds for i in 1:12
            _Fl_buf[i] *= dofs[i]
        end

        Pstart  = -_Fl_buf[1]
        Vystart =  _Fl_buf[2]
        Mystart =  _Fl_buf[6]
        Vzstart =  _Fl_buf[3]
        Mzstart = -_Fl_buf[5]

        fill!(P_tmp, Pstart)
        fill!(Vy_tmp, Vystart)
        fill!(Vz_tmp, Vzstart)
        @inbounds for j in 1:per_elem_res
            My_tmp[j] = Vystart * xinc[j] - Mystart
            Mz_tmp[j] = Vzstart * xinc[j] - Mzstart
        end

        for load in loadids
            accumulate_force!(load, xinc, P_tmp, My_tmp, Vy_tmp, Mz_tmp, Vz_tmp, ep)
        end

        x_offset = k > 0 ? xstore[k] : 0.0
        @inbounds for j in 1:per_elem_res
            k += 1
            xstore[k]  = x_offset + xinc[j]
            pstore[k]  = P_tmp[j]
            mystore[k] = My_tmp[j]
            vystore[k] = Vy_tmp[j]
            mzstore[k] = Mz_tmp[j]
            vzstore[k] = Vz_tmp[j]
        end
    end

    resize!(xstore, k); resize!(pstore, k); resize!(mystore, k)
    resize!(vystore, k); resize!(mzstore, k); resize!(vzstore, k)

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
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)
    
    element_loads = get_elemental_loads(model)

    IF = [ElementInternalForces(e, model.u[e.globalID], element_loads[e.elementID], inc_m, etype2DOF[typeof(e)]) for e in model.elements]

    # Preallocate concatenated arrays (avoids vcat + splat allocations)
    total_pts = sum(length(IF[i].x) for i in eachindex(IF))
    L_total = sum(IF[i].L for i in eachindex(IF))
    x  = Vector{Float64}(undef, total_pts)
    P  = Vector{Float64}(undef, total_pts)
    My = Vector{Float64}(undef, total_pts)
    Vy = Vector{Float64}(undef, total_pts)
    Mz = Vector{Float64}(undef, total_pts)
    Vz = Vector{Float64}(undef, total_pts)

    offset = 0
    x_cum = 0.0
    @inbounds for i in eachindex(IF)
        eif = IF[i]
        n = length(eif.x)
        x_off = x_cum + eif.L / 2
        for j in 1:n
            k = offset + j
            x[k]  = eif.x[j] + x_off
            P[k]  = eif.P[j]
            My[k] = eif.My[j]
            Vy[k] = eif.Vy[j]
            Mz[k] = eif.Mz[j]
            Vz[k] = eif.Vz[j]
        end
        offset += n
        x_cum += eif.L
    end

    return ElementInternalForces(nothing, L_total, x, P, My, Vy, Mz, Vz)
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
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)

    if isempty(model.frame_elements)
        return ElementInternalForces[]
    end

    results = Vector{ElementInternalForces}()
    sizehint!(results, length(model.frame_elements))
    element_loads_map = get_elemental_loads(model)

    for element in model.frame_elements
        L = ustrip(u"m", element.length)
        n = max(Int(round(L / inc_m)), 2)

        push!(results, ElementInternalForces(element, element_loads_map[element.elementID]; resolution = n))
    end

    return results
end


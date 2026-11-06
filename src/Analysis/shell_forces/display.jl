# ============================================================================
# Display
# ============================================================================

function Base.show(io::IO, sif::ShellInternalForces)
    print(io, "ShellInternalForces(")
    print(io, "N=[$(_fmt(sif.Nxx, 1)), $(_fmt(sif.Nyy, 1)), $(_fmt(sif.Nxy, 1))], ")
    print(io, "M=[$(_fmt(sif.Mxx, 2)), $(_fmt(sif.Myy, 2)), $(_fmt(sif.Mxy, 2))], ")
    print(io, "Q=[$(_fmt(sif.Qxz, 1)), $(_fmt(sif.Qyz, 1))])")
end

function Base.show(io::IO, ::MIME"text/plain", sif::ShellInternalForces)
    id = hasproperty(sif.element, :id) ? sif.element.id : "unknown"
    println(io, "ShellInternalForces for element ':$(id)'")
    println(io, "├─ Membrane forces [N/m]:")
    println(io, "│    Nxx = $(_fmt(sif.Nxx, 2))")
    println(io, "│    Nyy = $(_fmt(sif.Nyy, 2))")
    println(io, "│    Nxy = $(_fmt(sif.Nxy, 2))")
    println(io, "├─ Bending moments [N·m/m]:")
    println(io, "│    Mxx = $(_fmt(sif.Mxx, 4))")
    println(io, "│    Myy = $(_fmt(sif.Myy, 4))")
    println(io, "│    Mxy = $(_fmt(sif.Mxy, 4))")
    println(io, "└─ Transverse shear [N/m]:")
    println(io, "     Qxz = $(_fmt(sif.Qxz, 2))")
    print(io, "     Qyz = $(_fmt(sif.Qyz, 2))")
end

function Base.show(io::IO, sd::ShellDisplacements)
    id = hasproperty(sd.element, :id) ? sd.element.id : "unknown"
    print(io, "ShellDisplacements(:$id, w=$(_fmt(sd.w*1000, 3))mm)")
end

function Base.show(io::IO, ::MIME"text/plain", sd::ShellDisplacements)
    id = hasproperty(sd.element, :id) ? sd.element.id : "unknown"
    println(io, "ShellDisplacements for element ':$(id)'")
    println(io, "├─ Translations [m]:")
    println(io, "│    u = $(_fmt(sd.u, 6))")
    println(io, "│    v = $(_fmt(sd.v, 6))")
    println(io, "│    w = $(_fmt(sd.w, 6))")
    println(io, "└─ Rotations [rad]:")
    println(io, "     θx = $(_fmt(sd.θx, 6))")
    println(io, "     θy = $(_fmt(sd.θy, 6))")
    print(io, "     θz = $(_fmt(sd.θz, 6))")
end


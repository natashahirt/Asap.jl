# =============================================================================
# Element Processing
# =============================================================================

"""
    process_elements!(model::FrameModel)
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a model.
"""
function process_elements!(model::FrameModel)
    for element in model.elements
        fill!(element.Q, 0.0)
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    process_elements!(model::ShellModel)
    
Process shell elements in a shell model.
Thread-parallel: each element's LCS/R/K is independent.
"""
function process_elements!(model::ShellModel)
    elems = model.elements
    Threads.@threads for i in eachindex(elems)
        process!(elems[i])
    end
end

"""
    process_elements!(model::Model)
    
Process all elements in a unified model (both frame and shell).
Shell element processing is thread-parallel (each element is independent).
"""
function process_elements!(model::Model)
    # Process frame elements (typically few — sequential)
    for element in model.frame_elements
        fill!(element.Q, 0.0)
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
    
    # Process shell elements — embarrassingly parallel
    shells = model.shell_elements
    Threads.@threads for i in eachindex(shells)
        process!(shells[i])
    end
end

"""
    process_elements!(elements::Vector{<:FrameElement})
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a vector of elements.
"""
function process_elements!(elements::Vector{T}) where {T<:FrameElement}
    for element in elements
        fill!(element.Q, 0.0)
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    process_elements!(model::TrussModel)
    
Populate the transformation matrix and global elemental stiffness matrix of the elements in a truss model.
"""
function process_elements!(model::TrussModel)
    for element in model.elements
        lcs!(element, element.Ψ)
        R!(element)
        length!(element)
        global_K!(element)
    end
end

"""
    process_elements!(elements::Vector{<:ShellElement})

Process a vector of shell elements: compute LCS, transformation R, stiffness K, and mass M.
Thread-parallel: each element's processing is independent.
"""
function process_elements!(elements::Vector{T}) where {T<:ShellElement}
    Threads.@threads for i in eachindex(elements)
        process!(elements[i])
        populate_globalID!(elements[i])
    end
end


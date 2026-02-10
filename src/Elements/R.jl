"""
    _build_R12(Λ::NTuple{9,Float64}) -> Matrix{Float64}

Build the [12 × 12] block-diagonal rotation matrix from a flat 9-tuple
representing the 3×3 direction cosine matrix Λ (row-major).
Avoids 8 `zeros(3,…)` allocations from the old `[Λ zeros(3,9); …]` pattern.
"""
@inline function _build_R12(Λ::NTuple{9,Float64})
    R = zeros(12, 12)
    @inbounds for blk in (0, 3, 6, 9)          # four 3×3 blocks on the diagonal
        R[blk+1, blk+1] = Λ[1]; R[blk+1, blk+2] = Λ[2]; R[blk+1, blk+3] = Λ[3]
        R[blk+2, blk+1] = Λ[4]; R[blk+2, blk+2] = Λ[5]; R[blk+2, blk+3] = Λ[6]
        R[blk+3, blk+1] = Λ[7]; R[blk+3, blk+2] = Λ[8]; R[blk+3, blk+3] = Λ[9]
    end
    return R
end

"""
    R(element::Element; tol = 1e-4)

Get the [12 × 12] transformation matrix for a given element.

`tol` defines the threshold criteria for triggering special `R` computation for elements aligned with the global Y axis
"""
function R(element::Element; tol = 1e-4)
    
    xvec = local_x(element)
    CXx, CYx, CZx = xvec # local x cosines

    cΨ = cos(element.Ψ)
    sΨ = sin(element.Ψ)


    if norm(cross(xvec, globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = (0., CYx, 0.,
            -CYx*cΨ, 0., sΨ,
            CYx*sΨ, 0., cΨ)
    else # all other
        b1 = (-CXx * CYx * cΨ - CZx * sΨ) / sqrt(CXx^2 + CZx^2)
        b2 = sqrt(CXx^2 + CZx^2) * cΨ
        b3 = (-CYx * CZx * cΨ + CXx * sΨ) / sqrt(CXx^2 + CZx^2)

        c1 = (CXx * CYx * sΨ - CZx * cΨ) / sqrt(CXx^2 + CZx^2)
        c2 = -sqrt(CXx^2 + CZx^2) * sΨ
        c3 = (CYx * CZx * sΨ + CXx * cΨ) / sqrt(CXx^2 + CZx^2)

        Λ = (CXx, CYx, CZx, 
            b1, b2, b3, 
            c1, c2, c3)
    end
    
    return _build_R12(Λ)
end

function R!(element::Element; tol = 1e-6)

    CXx, CYx, CZx = element.LCS[1]
    cΨ = cos(element.Ψ)
    sΨ = sin(element.Ψ)

    if norm(cross(element.LCS[1], globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = [0. CYx 0.;
            -CYx*cΨ 0 sΨ;
            CYx*sΨ 0 cΨ]
    else # all other
        b1 = (-CXx * CYx * cΨ - CZx * sΨ) / sqrt(CXx^2 + CZx^2)
        b2 = sqrt(CXx^2 + CZx^2) * cΨ
        b3 = (-CYx * CZx * cΨ + CXx * sΨ) / sqrt(CXx^2 + CZx^2)

        c1 = (CXx * CYx * sΨ - CZx * cΨ) / sqrt(CXx^2 + CZx^2)
        c2 = -sqrt(CXx^2 + CZx^2) * sΨ
        c3 = (CYx * CZx * sΨ + CXx * cΨ) / sqrt(CXx^2 + CZx^2)

        Λ = [CXx CYx CZx; 
            b1 b2 b3; 
            c1 c2 c3]
    end

    @views element.R[1:3, 1:3] = Λ
    @views element.R[4:6, 4:6] = Λ
    @views element.R[7:9, 7:9] = Λ
    @views element.R[10:12, 10:12] = Λ

end

"""
    R(xvec::Vector{Float64}, Ψ; tol = 1e-4)

Return the [12 × 12] transformation matrix given a local x vector `xvec`

`tol` defines the threshold criteria for triggering special `R` computation for elements aligned with the global Y axis
"""
function R(xvec::Vector{Float64}, Ψ; tol = 1e-4)

    CXx, CYx, CZx = normalize(xvec) # local x cosines

    cΨ = cos(Ψ)
    sΨ = sin(Ψ)


    if norm(cross(xvec, globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = (0., CYx, 0.,
            -CYx*cΨ, 0., sΨ,
            CYx*sΨ, 0., cΨ)
    else # all other
        b1 = (-CXx * CYx * cΨ - CZx * sΨ) / sqrt(CXx^2 + CZx^2)
        b2 = sqrt(CXx^2 + CZx^2) * cΨ
        b3 = (-CYx * CZx * cΨ + CXx * sΨ) / sqrt(CXx^2 + CZx^2)

        c1 = (CXx * CYx * sΨ - CZx * cΨ) / sqrt(CXx^2 + CZx^2)
        c2 = -sqrt(CXx^2 + CZx^2) * sΨ
        c3 = (CYx * CZx * sΨ + CXx * cΨ) / sqrt(CXx^2 + CZx^2)

        Λ = (CXx, CYx, CZx, 
            b1, b2, b3, 
            c1, c2, c3)
    end
    
    return _build_R12(Λ)
end

"""
    R(element::TrussElement)

Get the [2 × 6] transformation matrix for a truss element.
"""
function R(element::TrussElement)
    Cx, Cy, Cz = element.LCS[1]
    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end

function R!(element::TrussElement)
    element.R[1, 1:3] = element.LCS[1]
    element.R[2, 4:6] = element.LCS[1]
end

"""
    R(xvec::Vector{Float64})

Get the [2 × 6] transformation matrix given a local x vector `xvec`
"""
function R(xvec::Vector{Float64})
    Cx, Cy, Cz = normalize(xvec)
    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end

function R(vec::SubArray)
    Cx, Cy, Cz = normalize(vec)
    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end


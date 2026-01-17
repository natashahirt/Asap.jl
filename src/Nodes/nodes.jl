abstract type AbstractNode end

"""
    Node(position::Vector{Quantity}, dofs::Vector{Bool}, id::Symbol = :node)
    Node(position::Vector{Quantity}, fixity::Symbol, id::Symbol = :node)

Instantiate a 6 DOF node with given position and fixities.

# Arguments
- `position`: Vector of 3 distance quantities [x, y, z]
- `dofs`: Vector of 6 booleans for DOF fixities, or a `fixity` symbol
- `id`: Optional symbol identifier

Available boundary conditions:
- :free, :fixed, :pinned
- :(x/y/z)free, :(x/y/z)fixed

# Example
```julia
Node([4.3u"m", 2.2u"m", 10.4u"m"], [true, true, false, true, false, false])
Node([4.3u"m", 2.2u"m", 10.4u"m"], :zfixed)
```
"""
mutable struct Node <: AbstractNode
    position::Vector{QuantityDistance}
    dof::Vector{Bool}
    nodeID::Int64
    globalID::Vector{Int64}
    reaction::Vector{Quantity}  # Forces (N) for DOFs 1-3, moments (N*m) for DOFs 4-6
    displacement::Vector{Quantity}  # Translations (m) for DOFs 1-3, rotations (rad) for DOFs 4-6
    id::Symbol

    function Node(position::Base.AbstractVector{<:Unitful.Quantity}, dofs::Vector{Bool}, id = :node)
        @assert length(position) == 3 && length(dofs) == 6 "Position vector must be in R³, DOFs must be length 6"
        
        # Convert all positions to meters
        pos_si = [uconvert(u"m", p) for p in position]

        node = new(
            pos_si,
            dofs,
            0,
            Vector{Int64}(undef, 6),
            Vector{Quantity}(undef, 6),  # Will be populated in postprocessing
            Vector{Quantity}(undef, 6),  # Will be populated in postprocessing
            id
        )

        return node
    end

    function Node(position::Base.AbstractVector{<:Unitful.Quantity}, fixity::Symbol, id = :node)
        @assert length(position) == 3 "Position vector must be in R³"

        dofs = copy(fixDict[fixity])
        
        # Convert all positions to meters
        pos_si = [uconvert(u"m", p) for p in position]

        node = new(
            pos_si,
            dofs,
            0,
            Vector{Int64}(undef, 6),
            Vector{Quantity}(undef, 6),  # Will be populated in postprocessing
            Vector{Quantity}(undef, 6),  # Will be populated in postprocessing
            id
        )

        return node
    end
end

"""
    TrussNode(position::Vector{Quantity}, dofs::Vector{Bool}, id::Symbol = :node)
    TrussNode(position::Vector{Quantity}, fixity::Symbol, id::Symbol = :node)

Instantiate a 3 DOF node with given position and fixities.

# Arguments
- `position`: Vector of 3 distance quantities [x, y, z]
- `dofs`: Vector of 3 booleans for DOF fixities, or a `fixity` symbol
- `id`: Optional symbol identifier

Available boundary conditions:
- :free, :pinned
- :(x/y/z)free, :(x/y/z)fixed

# Example
```julia
TrussNode([1.0u"m", 1.0u"m", 56.0u"m"], [false, true, true])
TrussNode([1.0u"m", 1.0u"m", 56.0u"m"], :pinned)
```
"""
mutable struct TrussNode <: AbstractNode
    position::Vector{QuantityDistance}
    dof::Vector{Bool}
    nodeID::Int64
    globalID::Vector{Int64}
    reaction::Vector{Quantity}  # All forces (N)
    displacement::Vector{Quantity}  # All translations (m)
    id::Symbol

    function TrussNode(position::Base.AbstractVector{<:Unitful.Quantity}, dofs::Vector{Bool}, id = :node)
        @assert length(position) == length(dofs) == 3  "Position and dof vector must be in R³"
        
        # Convert all positions to meters
        pos_si = [uconvert(u"m", p) for p in position]

        node = new(
            pos_si,
            dofs,
            0,
            Vector{Int64}(undef, 3),
            Vector{Quantity}(undef, 3),  # Will be populated in postprocessing (all forces in N)
            Vector{Quantity}(undef, 3),  # Will be populated in postprocessing (all translations in m)
            id
        )

        return node
    end

    function TrussNode(position::Base.AbstractVector{<:Unitful.Quantity}, fixity::Symbol, id = :node)
        @assert length(position) == 3 "Position vector must be in R³"

        dofs = copy(fixDict[fixity][1:3])
        
        # Convert all positions to meters
        pos_si = [uconvert(u"m", p) for p in position]

        node = new(
            pos_si,
            dofs,
            0,
            Vector{Int64}(undef, 3),
            Vector{Quantity}(undef, 3),  # Will be populated in postprocessing (all forces in N)
            Vector{Quantity}(undef, 3),  # Will be populated in postprocessing (all translations in m)
            id
        )

        return node
    end
end

"""
Common fixity types
"""
const fixDict = Dict(:fixed => [false, false, false, false, false, false],
    :free => [true, true, true, true, true, true],
    :xfixed => [false, true, true, true, true, true],
    :yfixed => [true, false, true, true, true, true],
    :zfixed => [true, true, false, true, true, true],
    :xfree => [true, false, false, false, false, false],
    :yfree => [false, true, false, false, false, false],
    :zfree => [false, false, true, false, false, false],
    :pinned => [false, false, false, true, true, true])

"""
Inactive DOF w/r/t plane
"""
const planeDict = Dict(:XY => [3, 4, 5],
    :YZ => [1, 5, 6],
    :ZX => [2, 4, 6])


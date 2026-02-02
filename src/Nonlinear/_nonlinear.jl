#=
Nonlinear Analysis Module for Asap.jl
=====================================

Provides geometric nonlinearity, buckling, and nonlinear static analysis.

Analysis Types:
- solve_static!     : Linear static (alias for solve!)
- solve_pdelta!     : P-delta iteration with geometric stiffness
- solve_buckling!   : Linear buckling eigenvalue problem
- solve_nonlinear!  : Full Newton-Raphson nonlinear statics

References:
- Cook, Malkus, Plesha "Concepts and Applications of FEA"
- McGuire, Gallagher, Ziemian "Matrix Structural Analysis"
- FinEtoolsFlexStructures.jl by Petr Krysl (MIT License) - geometric stiffness formulations
=#

# Result types
include("types.jl")

# Geometric stiffness assembly
include("geometric_stiffness.jl")

# Solvers
include("buckling.jl")
include("pdelta.jl")
include("nonlinear_statics.jl")

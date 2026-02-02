module Asap

using LinearAlgebra, SparseArrays
using Unitful
using StructuralBase: StructuralUnits  # Shared unit definitions (lbf, kip, ksi, psf)

# Units module (for type aliases and compatibility helpers)
include("Units/units.jl")
include("Units/compat.jl")

# global axes
const globalX::Vector{Float64} = [1., 0., 0.]
const globalY::Vector{Float64} = [0., 1., 0.]
const globalZ::Vector{Float64} = [0., 0., 1.]

include("Materials_Sections/material.jl")
include("Materials_Sections/section.jl")
include("Materials_Sections/layup.jl")  # Composite laminate definitions
export Material
export ShellMaterial
export Concrete_Shell
export Steel_Shell
export Section
export TrussSection
export Steel_Nmm
export Steel_kNm
# Section helpers
export is_timoshenko
export is_bernoulli_euler
# Composite materials (core infrastructure only - define presets in your project)
export Ply
export Laminate
export thickness
export isotropic_ply
export laminate_stiffnesses
export laminate_transverse_shear_stiffness
export laminate_inertias
export symmetric_laminate

include("Nodes/nodes.jl")
include("Nodes/utilities.jl")
export Node
export TrussNode
export planarize!
export fixnode!

include("Elements/elements.jl")
include("Elements/K.jl")
include("Elements/R.jl")
include("Elements/shell/shell.jl")       # Full shell elements (membrane + bending)
include("Elements/shell/composite.jl")   # Composite laminate shells
include("Elements/shell/meshing.jl")     # Shell meshing utilities
include("Elements/spring.jl")            # Grounded spring elements
include("Elements/utilities.jl")
export Element
export BridgeElement
export TrussElement
export has_eccentricity
export ShellElement
export ShellTri3
export CompositeShellTri3
export Spring
export spring_stiffness
export translational_stiffness
export rotational_stiffness
# Shell creation
export ShellSection
export Concrete_ShellSection_150mm
export Concrete_ShellSection_200mm
export Shell
export mesh
export get_nodes
# Diaphragm (rigid in-plane, massless)
export DiaphragmSection
export RigidDiaphragm
export Diaphragm
export stress
export bending_moments
export membrane_forces
export ply_stresses
export release!
export endpoints
export midpoint
export axial_force

# DYNAMICS - Mass matrices for modal/dynamic analysis
include("Dynamics/mass.jl")
export MassMatrixType
export MASS_CONSISTENT
export MASS_CONSISTENT_NO_ROTATION
export MASS_LUMPED
export MASS_LUMPED_NO_ROTATION
export local_mass
export global_mass
export global_mass!

include("Loads/loads.jl")
include("Loads/utilities.jl")
include("Loads/fixed_end_forces.jl")
include("Loads/shell_loads.jl")
export AbstractLoad
export NodeForce
export NodeMoment
export LineLoad
export GravityLoad
export PointLoad
export TributaryLoad
export AreaLoad
export SelfWeight
export SurfaceLoad        # Alias for backward compatibility
export intensities
export nodal_forces

include("Model/model.jl")
include("Model/utilities.jl")
# Model types
export AbstractModel
export ElementModel    # Base for single-element-type models
export Model           # Unified model (frames + shells)
export FrameModel      # Frame elements only (FrameModel <: ElementModel)
export ShellModel      # Shell elements only (ShellModel <: ElementModel)
export TrussModel      # Truss elements only (TrussModel <: ElementModel)
# Model helpers
export has_frame_elements
export has_shell_elements
export is_mixed
export all_elements
export n_elements
export update_DOF!
export connectivity
export node_positions
export volume


include("Model/preprocessing.jl")
include("Model/bridgeprocessing.jl")
include("Model/postprocessing.jl")
include("Model/analysis.jl")
export process!
export solve!
export solve
export to_displacement_vec
export to_reaction_vec
export assemble_stiffness
export assemble_stiffness!
export populate_globalID!
export add_springs!

# DYNAMICS - Modal analysis (requires Model to be defined)
include("Dynamics/modal.jl")
export ModalResult
export modal_analysis
export modal!
export modal
export natural_frequencies
export mode_shapes
export assemble_mass_matrix
export assemble_mass_matrix!
export print_modal_summary

# FORCE DENSITY METHOD
include("FDM/FDM.jl")
export FDMnode
export FDMelement
export FDMload
export Network
export force
export vector
export forces
export initial_lengths
export update_q!

# ANALYSIS - Internal forces and displacements
include("Analysis/translations.jl")
include("Analysis/force_functions.jl")
include("Analysis/force_analysis.jl")
include("Analysis/displacements.jl")
include("Analysis/shell_forces.jl")
export etype2DOF
export planarDOFs
export groupbyid
export ElementInternalForces
export ForceEnvelopes
export load_envelopes
export get_elemental_loads
export ElementDisplacements
export displacements
# Shell internal forces
export ShellInternalForces
export principal_moments
export principal_forces
export von_mises_stress
export max_surface_stresses
# Shell displacements
export ShellDisplacements
export max_deflection
# Unified interface (dispatches to ElementInternalForces or ShellInternalForces)
export InternalForces
export Displacements

# TRIBUTARY AREA COMPUTATION
include("Tributary/_tributary.jl")
# Types
export TributaryPolygon
export TributaryBuffers
export VertexTributary
export SpanInfo
# Edge tributaries (straight skeleton / one-way)
export get_tributary_polygons
export get_tributary_polygons_isotropic
export get_tributary_polygons_one_way
export vertices  # for converting parametric → absolute coords
# Vertex tributaries (Voronoi)
export compute_voronoi_tributaries
# Span calculations
export get_polygon_span
export governing_spans
export short_span, long_span, two_way_span

end 

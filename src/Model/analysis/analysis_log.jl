# =============================================================================
# Analysis log — opt-in factorization / solve accounting
# =============================================================================
#
# Counts every factorization and RHS solve, tagged by a caller-set stage label
# and the model it ran on, so a profiler can split a whole-building design run
# by stage and by model class (building frame, FEA slab mesh, story shell, …).
#
# The whole thing hangs off one `ScopedValue`. When no log is installed, every
# hook is a single `ScopedValue` read returning `nothing`, so production cost is
# nil and no function signature changes. `ScopedValue`s are inherited by spawned
# tasks, so labels stay correct across the parallel slab / MIP loops; the lock
# guards the shared event vector those tasks push into.

using Base.ScopedValues: ScopedValue, with

"""
    AnalysisEvent

One factorization, assembly or solve. `label` is the caller's stage path (e.g.
`"members/fp_column_check"`), `model_id` is `objectid(model)`, `kind` is one of
`:assemble`, `:factor`, `:refactor`, `:rhs`, `:dense_factor`, `n_free` is the
size of the factored system and `n_rhs` the number of right-hand sides solved.
"""
struct AnalysisEvent
    label::String
    model_id::UInt
    kind::Symbol
    n_free::Int
    n_rhs::Int
end

"""
    AnalysisLog

Thread-safe collector of [`AnalysisEvent`](@ref)s. Install one with
[`with_analysis_log`](@ref) (or by setting [`ANALYSIS_LOG`](@ref) directly) to
turn accounting on for the dynamic extent of a call.
"""
mutable struct AnalysisLog
    lock::ReentrantLock
    events::Vector{AnalysisEvent}
end
AnalysisLog() = AnalysisLog(ReentrantLock(), AnalysisEvent[])

"The installed [`AnalysisLog`](@ref), or `nothing` when accounting is off."
const ANALYSIS_LOG = ScopedValue{Union{Nothing,AnalysisLog}}(nothing)

"The current stage label applied to new [`AnalysisEvent`](@ref)s."
const ANALYSIS_LABEL = ScopedValue{String}("unlabeled")

"""
    log_analysis!(model, kind; n_free = 0, n_rhs = 0)

Record one [`AnalysisEvent`](@ref) under the current [`ANALYSIS_LABEL`](@ref).
No-op (one `ScopedValue` read) when no [`AnalysisLog`](@ref) is installed.
"""
@inline function log_analysis!(model, kind::Symbol; n_free::Int = 0, n_rhs::Int = 0)
    log = ANALYSIS_LOG[]
    log === nothing && return nothing
    ev = AnalysisEvent(ANALYSIS_LABEL[], objectid(model), kind, n_free, n_rhs)
    lock(log.lock) do
        push!(log.events, ev)
    end
    return nothing
end

"""
    with_analysis_log(f, log::AnalysisLog)

Run `f` with `log` installed as the active [`ANALYSIS_LOG`](@ref).
"""
with_analysis_log(f, log::AnalysisLog) = with(f, ANALYSIS_LOG => log)

"""
    with_analysis_label(f, label::AbstractString)

Run `f` with [`ANALYSIS_LABEL`](@ref) set to `label`.
"""
with_analysis_label(f, label::AbstractString) = with(f, ANALYSIS_LABEL => String(label))

# =============================================================================
# Per-thread scratch buffers
# =============================================================================
#
# A pool filled at the top level of a source file is built while the package is
# precompiled, which runs on one thread, and every later load reuses that
# array. Indexing it with `Threads.threadid()` then runs past the end. So pools
# are declared empty and filled by `_init_thread_pools!` from `Asap.__init__`.
# They are sized to `Threads.maxthreadid()`, not `nthreads()`: a process started
# with `-t auto` reports 8 threads here and hands out ids up to 16.

"""
    _fill_thread_pool!(pool, make) -> pool

Replace the contents of `pool` with one `make()` result per possible thread id.
"""
function _fill_thread_pool!(pool::Vector, make)
    empty!(pool)
    for _ in 1:Threads.maxthreadid()
        push!(pool, make())
    end
    return pool
end

"""
    _thread_buffer(pool)

The entry of `pool` that belongs to the calling thread.
"""
@inline _thread_buffer(pool::Vector) = pool[Threads.threadid()]

"""
    _init_thread_pools!()

Fill every per-thread pool in Asap. Called from `Asap.__init__`; a new pool must
be added here or it stays empty and the first `_thread_buffer` call throws.
"""
function _init_thread_pools!()
    _fill_thread_pool!(_K_WS_POOL, _KWorkspace)
    _fill_thread_pool!(_MOMENT_WS, ShellMomentWorkspace)
    _fill_thread_pool!(_SHELL_FORCES_WS, ShellForcesWorkspace)
    _fill_thread_pool!(_POPULATE_LOAD_BUF, () -> zeros(12))
    return nothing
end

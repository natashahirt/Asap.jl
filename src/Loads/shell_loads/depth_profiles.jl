

# =============================================================================
# Helper: Extract Depth Profile from Tributary Polygon
# =============================================================================

"""Internal: the tributary's load profile as `(positions, depths)` on `[0, 1]`."""
function _extract_depth_profile(tp)
    return copy(tp.profile_s), copy(tp.profile_w)
end

"""Record a skeleton arc from p_old to p_new separating face_left and face_right."""
function _record_skeleton_arc!(dcel::DCEL, registry::VertexRegistry, 
                                p_old::NTuple{2,Float64}, p_new::NTuple{2,Float64},
                                face_left::Int, face_right::Int)
    v1 = get_or_create_vertex!(registry, p_old)
    v2 = get_or_create_vertex!(registry, p_new)
    
    # Don't create self-loops
    v1 == v2 && return
    
    # Create halfedge pair marked as SKELETON edges (face separators):
    # Walking from v1 to v2, face_left is on the LEFT
    # DCEL convention: halfedge's face = the face on its LEFT
    # h1: v1 → v2 with face = face_left
    # h2: v2 → v1 with face = face_right (original right is now on left)
    # These edges SEPARATE faces and should not be crossed during face walks
    create_halfedge_pair!(dcel, v1, v2, face_left, face_right; is_skeleton=true)
end


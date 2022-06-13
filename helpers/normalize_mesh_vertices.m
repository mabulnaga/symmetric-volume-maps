function [mesh, scale] = normalize_mesh_vertices(mesh)
%NORMALIZE_MESH_VERTICES
%This function normalizes the mesh vertices and landmarks so that the
%resulting mesh has volume 1.
%Input: Mesh. struct with fields 'verts', 'tets', 'landmarks'.
%Output: Mesh with normalized vertices and landmarks.
%       scale: what scale was used for normalization.
V1 = sum(abs(volume(mesh.verts,mesh.tets)));
scale = (V1^(1/3));
mesh.verts = mesh.verts./scale;
mesh.orig_vol = V1;
mesh.normalization_scale = scale;
end


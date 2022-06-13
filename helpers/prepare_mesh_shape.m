function [Mesh] = prepare_mesh_shape(elem, node, surfaceOnly)
%PREPARE_MESH Function takes as input a 3D connectivity list and nodes and
% prepares the proper format for mapping.
%Inputs:
%       elem: NTx4: tetrahedron connectivity
%       node: Nx3: vertices in R^3
%       surfaceOnly: 1 if surface only, 0 otherwise.
T = elem;
V = node;
T = preprocess_flip_volume(T, V);
TR = triangulation(T, V);
E = edges(TR);
if(surfaceOnly == 1)
    [T, V] = freeBoundary(TR);
    face = T;
    E = edges(TR);
end
Mesh.verts = V;
Mesh.tets = T;
Mesh.edges = E;
Mesh.nv = length(V);
Mesh.nt = length(T);
Mesh.ne = length(E);
Fs = freeBoundary(TR);
if(surfaceOnly == 1)
    Mesh.boundary_faces = face;
else
    Mesh.boundary_faces = Fs;
end
end


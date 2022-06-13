function [Mesh] = prepare_mesh(img, create_mesh_param, surfaceOnly, varargin)
%PREPARE_MESH Function takes as input a 3D volume from an image and
% generates a volumetric or surface mesh, and a struct containing the
% proper format for mapping.
%Inputs:
%       img: 3D volume
%       create_mesh_params: parameters for generating the mesh (1x2 array)
%       surfaceOnly: 1 if surface only, 0 otherwise.
%       varargin: 1 to prepare a sphere mesh
%           2: rr (1x3)
%           3: elem (Ntx5)
%           4: node (Nvx3). If passing these in.
useSphere = 0;
if(nargin > 3)
    useSphere = varargin{1};
    rr = varargin{2};
end
if(useSphere == 0)
    [node,elem,face] = create_mesh(img,create_mesh_param(1), create_mesh_param(2));
else
    [node, face, elem] = meshanellip([0 0 0 ],rr,create_mesh_param(1), create_mesh_param(2));
end
T = elem(:,1:4);
V = node;
T = preprocess_flipVolume(T, V);
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


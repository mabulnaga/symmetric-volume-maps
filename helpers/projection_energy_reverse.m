function [E, grad, targetVerts] = projection_energy_reverse(X, Mesh1, Mesh2,weights, targetVerts)
% prepares and compute the reverse projection energy. This is from the
% target to the source. Usage: if the problem is map form mesh 1 to mesh 2,
% then call projection_energy_reverse(Mesh2.verts, Mesh2, Mesh1, Mesh2.A_1,
% []). Returns the gradient with respect to Mesh1.
% Inputs:
%       X: set of vertices to project. Assumes these are on mesh1.
%       Mesh1: source mesh
%       Mesh2: target mesh
%       weights: size Nvx1 of weights
%       targetVerts: optional. Location of vertices to map to.

if nargin == 4
    targetVerts = [];
end
% extract the boundary vertices
bdIndices = unique(Mesh1.boundary_faces(:));
Xb = X(bdIndices,:);
if(isempty(targetVerts))
    Xt = Mesh2.verts(unique(Mesh2.boundary_faces),:);
    [targetVerts, which_tri, bary_tri] = find_proj_target(Mesh2.boundary_faces,Xb,Xt);
end
[E,grad_proj] = compute_proj_energy_reverse(Xb,targetVerts,weights, Mesh2.boundary_faces, which_tri, bary_tri);
% now, reorganize the gradient to be with respect to all vertices.
grad = zeros(size(Mesh1.verts));
grad(bdIndices,:) = grad_proj;
end
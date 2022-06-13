function [E, grad, targetVerts] = projection_energy(X, Mesh1,Mesh2,gamma, weights)
% prepares and compute the projection energy.
% Inputs:
%       X: set of vertices to project. Assumes these are on mesh1.
%       Mesh1: source mesh
%       Mesh 2: target mesh
%       gamma: weight to use on the energy
%       weights: size Nvx1 of weights
%       targetVerts: optional. Location of vertices to map to. Removed. Now
%       it computes it each time.

% extract the boundary vertices
bdIndices = Mesh1.bd_indices;
Xb = X(bdIndices,:);
Xt = Mesh2.verts;
targetVerts = find_proj_target(Mesh2.boundary_faces,Xb,Xt);
%normalize the weights
w = weights(bdIndices)/Mesh1.s;
[E,grad_proj] = compute_proj_energy(Xb,targetVerts, w,gamma);
% now, reorganize the gradient to be with respect to all vertices.
grad = zeros(size(Mesh1.verts));
grad(bdIndices,:) = grad_proj;
end
function [E,grad] = projection_energy_reverse(X, Mesh1,Mesh2, gamma, weights)
% prepares and compute the reverse projection energy (from Mesh 2 to Mesh
% 1).
% Inputs:
%       X: set of vertices to target. Assumes these are on mesh1. (ex:
%       X_12) size Nv_1x1
%       Mesh1: source mesh
%       Mesh 2: target mesh
%       gamma: weights to use
%       weights: size Nv_2x1 of weights (barycentric area on Mesh 2).

bdIndices = Mesh2.bd_indices;
Xb = Mesh2.verts(bdIndices,:);
T1 = Mesh1.boundary_faces;
if(~isempty(weights))
    w = weights(bdIndices)/Mesh2.s;
else
    warning('no weights passed, using uniform');
    w = ones(length(Xb),1);
end
[E, grad_surface,idx] = compute_proj_energy_reverse(Xb,X,T1,w,gamma);
grad = zeros(length(X),3);
%idx = find(Mesh1.bd_indices);
inds = false(length(X),1);
inds(idx) = 1;
grad(inds,:) = grad_surface(inds,:);

end


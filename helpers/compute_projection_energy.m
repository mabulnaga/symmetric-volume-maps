function [E, grad] = compute_projection_energy(V, Vtarg, weights)
%Function computes the projection energy and returns the gradient.
%takes as input the set of vertices in mesh 1 and the target location of
%points. Also takes in the weight for each vertex.
%Inputs:
%       V: Nvx3: number of boundary vertices by 3
%       Vtarg: Nvx3: number of boundary vertices by 3
%       weights: Nvx1: barycentric area weights
D = (V-Vtarg);
P = weights.^(1/2)*(V-Vtarg);
E = trace(P*P');
if nargout > 1
    grad = 2*weights.*D;
end
end


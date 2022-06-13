function [E, grad] = compute_proj_energy(X,Xt,weights, gamma)
%Computes the projection energy
%Inputs:
%       X: nvx3: source vertex locations
%       Xt: ntx3: target vertex locations
%       weights: weights per vertex (Nvx1)
%       gamma: weight to use on the energy
P = X-Xt;
%E = gamma*trace(weights.*(P*P'));
E= gamma*sum(weights.*P.^2,'all');
if nargout > 1
    grad = 2*gamma*weights.*P;
end
end


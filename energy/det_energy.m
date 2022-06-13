function [E, grad] = det_energy(X, Mesh1, weight)
%DET_ENERGY: function to compute energy (1-detJ)^2
%determinant of the Jacobian.
%Inputs:
%       X: Nx3 vector of coordinate values
%       Mesh1: mesh structs
%       weight: Ntx1 per tetrahedron weight 

if(size(weight,2)>size(weight,1))
    weight = weight';
end

D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
%compute the Jacobian matrices
[~, J1] = generate_J_tets(Mesh1.tets,X,Mesh1.Oinv, false);

[detX1, grad1] = determinant_per_vertex(J1, Mesh1.tets, X);
% compute the gradient per vertex
E = squeeze((1-detX1).^2);
E = sum(weight.*E);
grad1 = -2*(1-detX1).*grad1.*reshape(weight,[1,1,length(weight)]);
grad = compute_grad_vertex(grad1, D, Mesh1.Oinv, Mesh1.tets);


end
function [E, grad_f] = Dirichlet_energy(S)
% Function computes the Dirichlet energy and gradient with respect to the
% singular values.
%Input:
%       S: DxN matrix of singular values. D is 2 for surface and 3 for
%       volume. N is the number of tetrahedra.
%Output:
%       E: 1xN SDE per tetrahedron
%       grad_f: DxN: gradient w.r.t. singular values

E = sum((S).^2,1) ;
grad_f = 2*(S);
end


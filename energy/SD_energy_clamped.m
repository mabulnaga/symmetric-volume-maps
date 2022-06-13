function [E, grad_f] = SD_energy_clamped(S)
% Function computes the symmetric Dirichlet energy and gradient with respect to the
% singular values.
%Input:
%       S: DxN matrix of singular values. D is 2 for surface and 3 for
%       volume. N is the number of tetrahedra.
%Output:
%       E: 1xN SDE per tetrahedron
%       grad_f: DxN: gradient w.r.t. singular values
thresh = 0.5;
m = (2*thresh -2*thresh.^(-3));
id_pos = S>=thresh;
id_neg = S<thresh;
E = zeros(size(S));
E_pos = (S(id_pos)).^2 + (S(id_pos)).^(-2);
E_neg = m*S(id_neg) + thresh.^2 + thresh.^(-2) - m*thresh;
E(id_pos) = E_pos;
E(id_neg) = E_neg;
E = sum(E,1);

% gradients
grad_f_pos =2*(S(id_pos))-2*S(id_pos).^(-3);
grad_f_neg = m;

grad_f = zeros(size(S));
grad_f(id_pos) = grad_f_pos;
grad_f(id_neg) = grad_f_neg;
end


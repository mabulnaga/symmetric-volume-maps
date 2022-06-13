function [E, grad_f] = ARAP_energy_negative_linear(S, varargin)
% Function computes the ARAP energy and gradient with respect to the
% singular values.
%Input:
%       S: DxN matrix of singular values. D is 2 for surface and 3 for
%       volume. N is the number of tetrahedra or triangles.
%       varargin: weights. 1xN scalar to weight each SVD contribution. If 
%                    none passed, then assume equal weight.
%Output:
%       E: 1xN ARAP energy per tetrahedron
%       grad_f: DxN: gradient w.r.t. singular values
thresh = 0.1;
id_pos = S>=thresh;
id_neg = S<thresh;
E_pos = (S(id_pos)-1).^2;
E_neg = (1-thresh)*abs((S(id_neg)-1));
E = zeros(size(S));
E(id_pos) = E_pos;
E(id_neg) = E_neg;
E = sum(E,1);
%now, compute the gradient...
grad_f_pos = 2*(S(id_pos)-1);
grad_f_neg = -1*(1-thresh);
%grad_f_neg = -10*ones(size(S(id_neg)));
grad_f = zeros(size(S));
grad_f(id_pos) = grad_f_pos;
grad_f(id_neg) = grad_f_neg;

end
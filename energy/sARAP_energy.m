function [E,grad_f] = sARAP_energy(S, varargin)
% Function computes the sARAP energy and gradient with respect to the
% singular values.
%Input:
%       S: DxN matrix of singular values. D is 2 for surface and 3 for
%       volume. N is the number of tetrahedra or triangles.
%       varargin: weights. 1xN scalar to weight each SVD contribution. If 
%                    none passed, then assume equal weight.
%Output:
%       E: 1xN ARAP energy per tetrahedron
%       TODO: grad_f: DxN: gradient w.r.t. singular values

E = 1/2*sum((S-1).^2,1) + 1/2*(prod(S).*(sum(S.^(-1)-1.^2,1)));
%grad_f = 2*(S-1);
grad_f = [];
end
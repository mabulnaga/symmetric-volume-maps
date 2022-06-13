function [E, grad] = compute_f(J,f,alpha, varargin)
%Function computes the distortion energy
%Inputs:
%       J: Jacobian marix, 3x3xN
%       f: pointer to function that operates on a matrix. Input to f is a
%       set of singular values, 3xN. f outputs a vector 1xN, and the
%       gradient with respect to the singular values, 3xN.
%       alpha: scalar, what to weight the energy and gradients by.
%       varargin:
%           weights: weight for the contribution of each tetrahedron. By
%               default all 1. size is 1xN
%Outputs:
% E: scalar, energy
% grad: 3x3xN: gradient with respect to the Jacobian matrix.
N = size(J,3);
dim = size(J,1);
if(nargin == 3)
    weights = ones(1,N);
else
    weights = varargin{1};
end

% compute the SVD
[U,S,V] = compute_signed_SVD_batch(J);
% compute the distortion energy
[E,grad_f] = f(S);
E = alpha*E.*weights;
if nargout < 2
    return;
end
%compute the gradient
grad_f = alpha*weights.*grad_f;

% compute the gradient per vertex.
grad = zeros(dim, dim, N);
diags = reshape(bsxfun(@plus,(0:N-1)*9,(1:4:9).'),3,1,[]);
grad(diags) = grad_f;
grad = reshape(grad,size(J,1),size(J,2),[]);
grad = pagemtimes(U,grad);
grad = pagemtimes(grad,'none',V,'transpose');

end



function [ grad ] = gpu_compute_grad_ARAP( T,alpha, Xvol,Oinv, U, S, V)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
D= [-1 -1 -1;1 0 0;0 1 0;0 0 1];

Id = eye(size(S,1),size(S,2));
Id = repmat(Id,1,1,size(S,3));
E = S - Id;

% now, gradient. dE/d\sigma is just 2(\sigma - 1)
E_grad = 2*E;

%now, we can find dE/DJ, as U*diag(gradE)*V^T
Vt = pagefun(@transpose, V);
gradJ = pagefun(@mtimes, U, E_grad);
gradJ = pagefun(@mtimes, gradJ, Vt);
%now, grad dE/dX is (U*diag(gradE)*V^T)*(D*(Z*D)^-1)^T)
%gradJ = pagefun(@transpose, gradJ);
map_term = pagefun(@mtimes, D, Oinv);
map_term = pagefun(@transpose, map_term);
grad = pagefun(@mtimes, gradJ, map_term);
% multiply by the gradient
grad = bsxfun(@times, grad, reshape(Xvol,[1,1,size(grad,3)]));
%multiply by the regularization term penalty
%everything above here was the first term of the summation, now will add
%the derivative of the volume term.
%Now, convert to per vertex gradient.
T = reshape(T',1,[]);
gradMat = reshape(grad,3,[]);
z1 = accumarray(T',gradMat(1,:)')';
z2 = accumarray(T',gradMat(2,:)')';
z3 = accumarray(T',gradMat(3,:)')';
grad = alpha* [z1;z2;z3]';

end


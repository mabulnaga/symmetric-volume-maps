function [totalD, U, S, V] = gpu_compute_distortion_ARAP(X, alpha, Xvol,Oinv)
%Computes the ARAP energy, weighted by volume. Outputs the energy
%(scalar), and tensors containing, U, S, V of the SVD per tetrahedron.
D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
J = gpu_compute_J(X,Oinv);

% compute the SVD tensor
S = gpuArray(zeros(size(J)));
U = S;
V = S;
for i = 1 : size(J,3)
    [u,s,v] = svd(J(:,:,i));
    S(:,:,i) = s;
    U(:,:,i) = u;
    V(:,:,i) = v;
end

% compute ARAP energy
Id = eye(size(S,1),size(S,2));
Id = repmat(Id,1,1,size(S,3));
E = S - Id;
% ARAP is sum (sigma_i -1)^2, or the Frobenius norm square
Et = pagefun(@transpose, E);
E_mul = pagefun(@mtimes, E, Et);

%now, need to reshape, and weight by volume. (to get the trace)
diags = reshape(bsxfun(@plus,(0:size(J,3)-1)*9,(1:4:9).'),3,1,[]);
D = sum(E_mul(diags));
D = reshape(D,1,length(D));
scaledD = gpuArray(Xvol).*D;
%alpha is the weight on regularization 
totalD = alpha * sum(scaledD);

end


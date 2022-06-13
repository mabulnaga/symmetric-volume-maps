function [E,grad] = optimize_tet_inversion(X,V,f,Mesh,weights, nring)
% Used to optimize tet inversions with fmincon
% Inputs:
%       X: (3N)x1: vector of mesh coordinates
%       V: Nx3: matrix of original mesh coordinates
%       f: distortion function handle
%       Mesh: mesh data structure
%       weights: Ntx1 vector of weights
%       nring: nring neighborhood to use
% Returns:
%       E: scalar 
%       grad: (3N)x1: vector of gradients
use_GPU = false;
alpha = 1;
signed_svd = true;
D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
T = Mesh.tets;
Nt = size(T,1);
%compute the Jacobian matrices
Xmap = reshape(X,[],3);
[vertID, tetID, tetNeibID] = get_inverted_tet_verts(T,V,Xmap, nring);
tetIdx = zeros(1,size(T,1)); 
tetIdx(tetNeibID) = 1;
tetIdx = logical(tetIdx);
[~, J1] = generate_J_tets(Mesh.tets(tetIdx,:),Xmap,Mesh.Oinv(:,:,tetIdx), use_GPU);
[ED1, gradJ1] = compute_f(J1, f, signed_svd, alpha, weights(tetIdx));
E = sum(ED1);
gradJ = zeros(3,3,Nt);
gradJ(:,:,tetIdx) = gradJ1;
grad1 = compute_grad_vertex(gradJ, D, Mesh.Oinv, T);
grad = grad1(:);

end


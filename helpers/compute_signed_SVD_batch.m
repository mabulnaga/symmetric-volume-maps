function [U,d,V,S] = compute_signed_SVD_batch(F)
%This function computes the signed SVD according to the paper
%Invertible Finite Elements for Robust Simulation of Large Deformation, by
%Irving et al.
%Inputs:    
%       F: nxnxN Jacobian matrix. n=3 for tetrahedra, n=2 for triangle.
%Outputs:
%       U: nxnxN rotation matrix (detU = 1)
%       d: nx1xN singular value matrix
%       V: nxnxN rotation amtrix (detU = 1)
%       S: nxnxN diagonal matrix of signed singular values. only returned
%       if nargout = 4.
n = size(F,3);
dim = size(F,1);
if(~ismac)
    [U,S,V] = batchSVD3x3Eigen(F);
    diags = reshape(bsxfun(@plus,(0:n-1)*9,(1:4:9).'),3,1,[]);
    d = squeeze(S(diags));
else
   [U,S,V] =  batchop('svd', F, 3);
   if(isa(V,'double'))
        V = permute(V,[2 1 3]);
   end
   d = S;
end

% Assign the sign to the smallest singular value.
detV = det3x3(V);
ind = detV < 0;
iF = find(ind);
V(:,3,iF) = -1 * V(:,3,iF);
d(3,iF) = -1*d(3,iF);
% negate the singular value if detU < 0.
detU = det3x3(U);
ind = detU < 0;
iF = find(ind);
U(:,3,iF) = -1 * U(:,3,iF);
d(3,iF) = -1*d(3,iF);
if(nargout == 4)
    S(3,3,iF) = d(3,iF);
end
end


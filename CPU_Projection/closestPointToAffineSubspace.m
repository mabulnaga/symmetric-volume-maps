function [projection,barycentric] = closestPointToAffineSubspace(simplices, query)
% Projects points onto affine subspaces spanned by simplices.
%
% Inputs:
%    simplices (m x d x k):  m k-simplices in R^d
%    query (n x d):  n query points in R^d
%
% Outputs:
%    projection (m x n x d):  projected query points
%    barycentric (m x n x k):  barycentric coordinates
%
% Only works for k <= 4

n = size(query,1);
d = size(query,2);
m = size(simplices,1);
assert( size(simplices,2) == d );
k = size(simplices,3);

if k == 1 % boring case
    barycentric = ones([m,n,1]);
    projection = sum(reshape(barycentric,m,n,1,k).*reshape(simplices,m,1,d,k),4);
    return
end

% Vectors along edges of simplices
edgeVecs = zeros(m,d,k-1);
for i=2:k
    edgeVecs(:,:,i-1) = simplices(:,:,i) - simplices(:,:,1);
end

% Set up the matrix of the normal equations
matrices = zeros(m,k-1,k-1);
for i=1:(k-1)
    for j=1:(k-1)
        matrices(:,i,j) = sum(squeeze(edgeVecs(:,:,i).*edgeVecs(:,:,j)),2);
    end
end

% Invert everything (m x k-1 x k-1)
inverses = invertSymmetricMatrices(matrices);

% Projection matrix (V'*V)^{-1} * V' (m x k-1 x d)
% projectionMatrices(i,j,p) = sum_l inv(i,j,l)*edgeVecs(i,p,l)
projectionMatrices = sum(reshape(inverses,m,k-1,1,k-1).*reshape(edgeVecs,m,1,d,k-1),4);

% Compute barycentric coordinates

% Apply projection matrix
% edgeProjection(i,j,q) = sum_l projectionMatrices(i,q,l) * ( query(j,l) - simplices(i,l,1) )
edgeProjection = sum(reshape(projectionMatrices,m,1,k-1,d).*(reshape(query,1,n,1,d)- reshape(simplices(:,:,1),m,1,1,d)),4);

% Convert to barycentric coordinates (currently in coordinates of edge vectors)
barycentric = zeros(m,n,k);
barycentric(:,:,1) = 1;
for i=1:(k-1)
    barycentric(:,:,i+1) = edgeProjection(:,:,i);
    barycentric(:,:,1) = barycentric(:,:,1) - edgeProjection(:,:,i);
end

% Embed in R^d
% projection(i,j,p) = sum_l barycentric(i,j,l)*simplices(i,p,l)
projection = sum(reshape(barycentric,m,n,1,k).*reshape(simplices,m,1,d,k),4);
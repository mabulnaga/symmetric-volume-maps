function [projection,whichSimplex,barycentric] = closestSimplex(simplices, query)
% Projects points onto affine subspaces spanned by simplices.
%
% Inputs:
%    simplices (m x d x k):  m k-simplices in R^d
%    query (n x d):  n query points in R^d
%
% Outputs:
%    projection (n x d):  projected query points
%    whichSimplex (n x 1):  index of closest simplex
%    barycentric (n x k):  barycentric coordinates
%
% Only works for k <= 4

m = size(simplices,1);
d = size(simplices,2);
k = size(simplices,3);
n = size(query,1);

% Project onto all simplices
[allprojections,allbarycentric,sqdistances] = closestPointInSimplex(simplices,query);

% Find index with smallest distance
[~,whichSimplex] = min(sqdistances,[],1);

% Slice through tensors
projection = zeros(n,d);
r = repmat( (1:n)', d , 1);
c = reshape(repmat( 1:d, n, 1 ), [], 1);
projection( sub2ind([n,d],r,c) ) = allprojections( sub2ind([m,n,d],repmat(whichSimplex',d,1),r,c) );

barycentric = zeros(n,k);
r = repmat( (1:n)', k , 1);
c = reshape(repmat( 1:k, n, 1 ), [], 1);
barycentric( sub2ind([n,k],r,c) ) = allbarycentric( sub2ind([m,n,k],repmat(whichSimplex',k,1),r,c) );
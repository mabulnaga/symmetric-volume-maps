function [projection,barycentric,sqdistances] = closestPointInSimplex(simplices, query)
% Projects points onto affine subspaces spanned by simplices.
%
% Inputs:
%    simplices (m x d x k):  m k-simplices in R^d
%    query (n x d):  n query points in R^d
%
% Outputs:
%    projection (m x n x d):  projected query points
%    barycentric (m x n x k):  barycentric coordinates
%    sqdistances (m x n):  squared distance to projection
%
% Only works for k <= 4

m = size(simplices,1);
d = size(simplices,2);
k = size(simplices,3);
n = size(query,1);

simplexFaces = nchoose(1:k); % enumerate 2^k possible faces

barycentric = zeros(m,n,k);
sqdistances = ones(m,n)*Inf;

for i=1:length(simplexFaces)
    face = simplexFaces{i};
    if isempty(face), continue; end 
    
    % Project onto affine subspace spanned by this face
    faceVertices = simplices(:,:,face);
    [affineProjection,affineBarycentric] = closestPointToAffineSubspace(faceVertices, query);
    
    % We  want to throw away projections whose barycentrics aren't in [0,1]
    badProjections = sum(affineBarycentric < 0 | affineBarycentric > 1,3) > 0;
    affineBarycentric(repmat(badProjections,1,1,length(face))) = Inf;
    affineProjection(repmat(badProjections,1,1,d)) = Inf;
    
    % Distance squared from each query to each simplex
    curDistSquared = sum( (affineProjection - reshape(query,1,n,d)).^2 ,3);
    
    % If it's closer, keep it
    replace = (curDistSquared < sqdistances);
    sqdistances(replace) = curDistSquared(replace);

    % Annoying indexing to keep the barycentric coordinates of the projection
    barycentric(repmat(replace,1,1,k)) = 0;
    [row,col] = find(replace);
    for j=1:length(face)
        idx = sub2ind([m n k],row,col,face(j)*ones(size(row)));
        idx2 = sub2ind([m n k],row,col,j*ones(size(row)));
        barycentric(idx) = affineBarycentric(idx2);
    end
end

% I'm lazy and recompute the projected point from the barycentric coordinates here
projection = sum(reshape(barycentric,m,n,1,k).*reshape(simplices,m,1,d,k),4);
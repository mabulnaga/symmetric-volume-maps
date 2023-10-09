function [projection,index,barycentric] = approximateClosestPointInSimplex(simplices, query, nsamples)
% Projects points onto simplices, uses knn approximation.
%
% Inputs:
%    simplices (m x d x k):  m k-simplices in R^d
%    query (n x d):  n query points in R^d
%
% Outputs:
%    projection (n x d):  projected query points
%    index (n x 1):  Index of the approximate closest simplex
%    barycentric (n x k):  barycentric coordinates
%
% Only works for k <= 4

m = size(simplices,1);
d = size(simplices,2);
k = size(simplices,3);
n = size(query,1);

% Generate random barycentric coordinates
barycentric = rand(nsamples, k);
barycentric(:,1) = 0;
barycentric(:,k+1) = 1;
barycentric = sort(barycentric,2);
barycentric = barycentric(:,2:end) - barycentric(:,1:(end-1));
barycentric(1:k,1:k) = eye(k); % corners of tet

% Generate points inside the simplices
samples = zeros(m*nsamples,d);
for i=1:nsamples 
    p = zeros(m,d);
    for dim=1:k, p = p + barycentric(i,dim)*squeeze(simplices(:,:,dim)); end
    samples( ((i-1)*m + 1):(i*m) ,:) = p;
end

% Knn search
index = mod(knnsearch(samples,query)-1,m)+1; % multiple sample points per tet

% Compute barycentric projection onto the knn simplex as a refinement
simplexFaces = nchoose(1:k); % enumerate 2^k possible faces
barycentric = zeros(n,k);
indexedSimplices = simplices(index,:,:); % n x d x k
bestDist = Inf*ones(n,1);
projection = zeros(n,d);

for i=1:length(simplexFaces)
    face = simplexFaces{i};
    if isempty(face), continue; end 
    
    subsimplices = indexedSimplices(:,:,face);
    
    if length(face) == 1
        % Just a point
        p = squeeze(subsimplices);
        dists = sum( (p-query).^2, 2);
        
        replace = bestDist > dists;
        barycentric(replace,:) = 0;
        barycentric(replace,face) = 1;
        bestDist(replace) = dists(replace);
        projection(replace,:) = p(replace,:);
    else
        % A simplex
        kk = length(face);
        
        % Subtract base point from query points
        shiftedqueries = query - squeeze(subsimplices(:,:,1));
        
        % Vectors along edges of simplex
        edgevecs = zeros(n,d,kk-1);
        for u=2:kk, edgevecs(:,:,u-1) = subsimplices(:,:,u)-subsimplices(:,:,1); end
        
        % Gram matrix of edge vectors and its inverse
        matrices = zeros(n,kk-1,kk-1);
        for u=1:(kk-1),for v=1:(kk-1),for w=1:d
            matrices(:,u,v) = matrices(:,u,v) + edgevecs(:,w,u).*edgevecs(:,w,v);
        end,end,end
        invmatrices = invertSymmetricMatrices(matrices);
        
        % Right hand side for projection problem
        rhs = zeros(n,kk-1);
        for u=1:(kk-1),for v=1:d, rhs(:,u) = rhs(:,u) + edgevecs(:,v,u).*shiftedqueries(:,v); end,end
        
        coeff = zeros(n,kk-1);
        for u=1:(kk-1),for v=1:(kk-1), coeff(:,u) = coeff(:,u) + invmatrices(:,u,v).*rhs(:,v); end,end
    
        bary = [1-sum(coeff,2) , coeff];
        validbary = (sum(bary < 0,2) == 0);
        bary(~validbary,:) = Inf;
        
        points = zeros(n,d);
        for u=1:d,for v=1:kk, points(:,u) = points(:,u) + bary(:,v).*subsimplices(:,u,v); end,end
        
        dists = sum((points-query).^2,2);
        replace = bestDist > dists;
        
        barycentric(replace,:) = 0;
        barycentric(replace,face) = bary(replace,:);
        bestDist(replace) = dists(replace);
        projection(replace,:) = points(replace,:);
    end
end

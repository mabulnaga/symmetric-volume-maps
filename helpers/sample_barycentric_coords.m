function [barycentric] = sample_barycentric_coords(nsamples, simplex_dim)
%SAMPLE_BARYCENTRIC_COORDS: function outputs a set of barycentric
%coordinates.
%Inputs: nsamples: number of samples desired
%       simplex_dim: 4 for tetrahedron, 3 for triangle, etc...
%Outputs: barycentric: barycentric coordinates. Also assigns 1 to each
%vertex.
k = simplex_dim;
barycentric = rand(nsamples, k);
barycentric(:,1) = 0;
barycentric(:,k+1) = 1;
barycentric = sort(barycentric,2);
barycentric = barycentric(:,2:end) - barycentric(:,1:(end-1));
barycentric(1:k,1:k) = eye(k); % corners of tet
end


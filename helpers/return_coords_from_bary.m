function [samples] = return_coords_from_bary(T,V,barycentric_coords)
%RETURN_COORDS_FROM_BARY
%Outputs a set of points based on sampled barycentric coordinates.
%Inputs: 
%       T: mxk simplex, k is the simplex size (e.g. 4 for tetrahedron)
%       V: nxd, n points in R^d.
%       barycentric_coords: lxk set of barycentric coordinates, where l is
%       the number of coordinates per tet.
%Outputs:
%       samples: (mxl)xd sampled points
m = size(T,1);
k = size(T,2);
n = size(V,1);
d = size(V,2);
nsamples = size(barycentric_coords,1);
np = m*nsamples;

samples = zeros(np,d);
% generate simplices
simplices = generate_tets(T,V,false);
simplices = permute(simplices,[3 1 2]);
for i=1:nsamples 
    p = zeros(m,d);
    for dim=1:k, p = p + barycentric_coords(i,dim)*squeeze(simplices(:,:,dim)); end
    %samples( ((i-1)*m + 1):(i*m) ,:) = p;
    samples((i : nsamples : np),:) = p;
end

end


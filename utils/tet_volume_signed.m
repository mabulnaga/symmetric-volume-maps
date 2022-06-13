function [ vol ] = tet_volume_signed( T )
% computes the signed volume of one tetrahedron.
% Inputs:
%       T: 3x4 matrix, each column is a vertex in R^3.
% Outputs:
%       vol: scalar with the signed volume.
%Assumes the tetrahedron is 3x4, where each vertex is a column vector.
a = ones(4,1);
vol = 1/6*(det([T',a]));
end


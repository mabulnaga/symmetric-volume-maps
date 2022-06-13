function [vertID, tetID, tetNeibID] = get_inverted_tet_verts(T,X,Xmap, nring,boundaryIndices)
%Gets the indices of vertices in inverted tets.
%Inputs:
%       T: input Triangulation, Ntx4
%       X: original vertex locaitons Nx3
%       Xmap: mapped vertex locations
%       nring: n-ring neighborhood. 0 ring means just the vertices.
%       boundary Indices (optional): surface triangulation.
%Outputs:
%       vertID: indices of vertices in inverted tets or their neighbors:
%           vector
%       tetID: ID of inverted tets (logical Ntx1)
%       tetNeibID: indices of all tetrahedra containing vertices in vertID.

if nargin == 4
    boundaryIndices =[];
end
ind_invert = get_inverted_tet_ID(T,X,Xmap,boundaryIndices);
tetID = ind_invert;
% vertex IDs
tetNeibID = find(tetID);
vertID = unique(T(tetNeibID(:),:));

if(nring >= 1)
    A = adjacency_matrix(T); 
    Ann = compute_knn_graph(A,nring);
    [I,J] = find(Ann);
    idx = ismember(I,vertID);
    J = unique(J(idx));
    vertID = unique([vertID(:);J(:)]);
    tetNeibId2 = [];
    V = vertexAttachments(triangulation(T,X),vertID);
    for i = 1 : size(V,1); tetNeibId2 =[tetNeibId2,V{i}];end
    tetNeibID = unique(tetNeibId2);
end
end



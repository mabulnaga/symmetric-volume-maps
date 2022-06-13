function [X_12, X_21] = inverted_tet_repair_split(f,X_1, X_2, X_12, X_21, Mesh1, Mesh2,nring, options)
%INVERTED_TET_REPAIR_SPLIT: wrapper function to optimize inverted
%tetrahedra.
%Inputs:
%       f: distortion energy function handle 
%       X_1: start vertex positions of mesh 1 N1x3
%       X_2: start vertex positions of mesh 2 N2x3
%       X_12: image of mesh 1 map: N1x3
%       X_21: image of mesh 2 map: N2x3
%       Mesh 1: struct of mesh 1
%       Mesh 2: struct of mesh 2
%       nring: number of tet neighbors to optimize over
%       options: fmincon options
tic
[~,tID_old] = get_inverted_tet_verts(Mesh1.tets,Mesh1.verts,X_12, nring);
[~,tID2_old] = get_inverted_tet_verts(Mesh2.tets,Mesh2.verts,X_21, nring);

vol_weight1 = Mesh1.vol_normalized;
vol_weight2 = Mesh2.vol_normalized;
% vol_weight1  =ones(size(vol_weight1));
% vol_weight2 = ones(size(vol_weight2));
% vol_weight1(tID_old) = 100;
% vol_weight2(tID2_old) = 100;
for i = 1 : 1
X = fmincon(@(X)optimize_tet_inversion(X, X_1, f, Mesh1, vol_weight1,nring),X_12(:),[],[],[],[],-Inf,Inf,[],options);
X_12 = reshape(X,[],3);
X = fmincon(@(X)optimize_tet_inversion(X, X_2, f, Mesh2, vol_weight2,nring),X_21(:),[],[],[],[],-Inf,Inf,[],options);
X_21 = reshape(X,[],3);
toc
end
[~,tID_p] = get_inverted_tet_verts(Mesh1.tets,Mesh1.verts,X_12, nring);
[~,tID_p2] = get_inverted_tet_verts(Mesh2.tets,Mesh2.verts,X_21, nring);
fprintf('*************** # inverted tets prev: %d, %d, # inverted tets after: %d , %d\n ',length(find(tID_old)), length(find(tID2_old)), length(find(tID_p)), length(find(tID_p2)));
end


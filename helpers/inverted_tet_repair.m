function [X_12, X_21] = inverted_tet_repair(X_1, X_2, X_12, X_21, P_12, P_21, Mesh1, Mesh2, c_1, c_2, nring, energy_repair, options, post_convergence_repair)
%inverted_tet_repair: repairs inverted tetrahedra by optimizing the
%distortion energy on vertices in flipped tetrahedra.
%Inputs:
%   X_1: n1x3 Mesh 1 vertices
%   X_2: n2x3 Mesh 2 vertices
%   X_12: n1x3 image of vertices of the map from Mesh1 to Mesh2
%   X_21: n2x3 image of vertices of the map from Mesh2 to Mesh1
%   P_12: n1xn2: map from Mesh1 to Mesh2
%   P_21: n2xn1: map from Mesh2 to Mesh1
%   c_1: total volume of Mesh 1
%   c_2: total volume of Mesh2
%   nring: n-ring neighborhood of vertices in inverted tetrahedra to
%       optimize
%   energy_repair: function handle of distortion energy to repair
%   options: fmincon options
%   post_convergence_repair: flag to determine if we are in the
%       post-convergence repair stage.
% Outputs:
%   X_12, X_21: new set of mapped vertex positions after inverted
%   tetrahedron repair.

count = 0;
alpha = 1;
beta = 0; gamma = 0;
% by default, set to ignore all vertices
ignore_verts{1} = unique(Mesh1.tets(:));
ignore_verts{2} = unique(Mesh2.tets(:));
vol_weight1 = Mesh1.vol_normalized;
vol_weight2 = Mesh2.vol_normalized;

%get the vertices for the tets to uninvert. if post_convergence_repair, keep doing
%until we have 0 inversions, for a maximum of 5 LBFGS calls.
while(true)
    [repair_verts{1},flipped_tet_ID] = get_inverted_tet_verts(Mesh1.tets,Mesh1.verts,X_12, nring);
    [repair_verts{2},flipped_tet_ID2] = get_inverted_tet_verts(Mesh2.tets,Mesh2.verts,X_21, nring);
    if(count == 0)
        flipped_tet_ID_old = flipped_tet_ID; 
        flipped_tet_ID2_old = flipped_tet_ID2;
    end
    num_inversions = length(repair_verts{1})+length(repair_verts{2});
    if(num_inversions>0 && count < 5)
        [X_min] = fmincon(@(X)compute_f_all(Mesh1, Mesh2, X, P_12, P_21, Mesh1.V_1, Mesh2.V_1, X_1, X_2, energy_repair,...
                c_1, c_2,  alpha, beta, gamma, vol_weight1,vol_weight2, Mesh1.A_1, Mesh2.A_1, ignore_verts, repair_verts),[X_12(:);X_21(:)],[],[],[],[],-Inf,Inf,[],options);
        %update X_12, X_21.
        X_12 = X_min(1:Mesh1.nv*3,:);
        X_12 = reshape(X_12,[],3);
        X_21 = X_min(Mesh1.nv*3+1:end,:);
        X_21 = reshape(X_21,[],3);
        count = count + 1;
    else
        break;
    end
    if(~post_convergence_repair)
        break;
    end
end
[~,flipped_tet_ID_p] = get_inverted_tet_verts(Mesh1.tets,Mesh1.verts,X_12, nring);
[~,flipped_tet_ID_p2] = get_inverted_tet_verts(Mesh2.tets,Mesh2.verts,X_21, nring);

fprintf('*************** # inverted tets prev: %d, %d, # inverted tets after: %d , %d\n ',length(find(flipped_tet_ID_old)), length(find(flipped_tet_ID2_old)), length(find(flipped_tet_ID_p)), length(find(flipped_tet_ID_p2)));

end


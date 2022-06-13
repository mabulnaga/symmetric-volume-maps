% P is a full functional map from M1 to M2 (M1.nv x M2.nv).
% T is a matrix with 4 columns, first is a vector with indices of target 
% faces per vertex of M1, 3 last hold target barycentric coordinates per vertex
% Written by Danielle Ezuz and modified by Mazdak Abulnaga.
function P = precise_T_to_P_tet(T, tets, nv)

P = sparse(repmat((1:size(T,1))', 4,1), ...
    reshape(double(tets(T(:,1),:)),[],1), ...
    reshape(T(:,2:5),[],1), ...
    size(T,1), nv);
    
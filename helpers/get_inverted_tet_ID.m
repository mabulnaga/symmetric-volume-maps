function [ind_invert,onBD_id, notBD_id, num_inv, num_bd, num_not_bd] = get_inverted_tet_ID(T,X, Xmap, boundaryIndices)
%This function finds the indices of all tetrahedra.
%Can also return the indices of the inverted tetrahedra with nodes on the
%boundary, and nodes not on the boundary.
%Inputs:
%       T: tetrahedra
%       X: original mesh vertices
%       boundaryIndices: indices of the boundary triangles ntx3
%       Xmap: mapped vertex locations
%Outputs:
%       ind_invert: logical for where the inverted tetrahedra are
%       onBD_id: id of tetrahedra on the boundary
%       notBD_id: id of tetrahedra in the interior
% example usage:
% get_inverted_tet_ID(Mesh.tets,Mesh.verts,Ts, P_12*Mesh2.verts);
if(nargin == 3)
    Ts = [];
else
    Ts = boundaryIndices;
end
det_J_forward = compute_distortion_J(T, T, X, Xmap);
ind_invert = squeeze(det_J_forward) <= 0;

if(nargout > 1)
    onBD = 0;
    onBD_id = [];
    notBD_id = [];
    invert_ids = find(ind_invert);
    if(~isempty(Ts))
        for i = 1 : length(invert_ids)
            ind = invert_ids(i);
            if( find(ismember(T(ind,:),unique(Ts(:)))))
                onBD = onBD + 1;
                onBD_id = [onBD_id, ind];
            else
                notBD_id = [notBD_id, ind];
            end
        end
    end
    num_inv = length(find(invert_ids));
    num_bd = length(find(onBD_id));
    num_not_bd = num_inv - num_bd;
end
end


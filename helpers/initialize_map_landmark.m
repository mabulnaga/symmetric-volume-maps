function [P12_init, P21_init] = initialize_map_landmark(Mesh1, Mesh2, landmarks1, landmarks2)
%INITIALIZE_MAP_LANDMARK
% Initializes the maps P_12 and P_21 based on a set of input landmarks.
% Inputs:
%       Mesh1, Mesh2: structs containing mesh information
%       landmarks1, landarks2: indices into Mesh vertices with matching
%           landmarks
%       initialized_map: boolean, whether the map has already been
%           initialized by a surface map or not
%       fix_landmarks:
    X_1 = Mesh1.verts;
    P12_init = sparse((1:Mesh1.nv)', landmarks2(knnsearch(X_1(landmarks1,:), X_1)), ones(Mesh1.nv,1), Mesh1.nv, Mesh2.nv);
    P21_init = sparse((1:Mesh2.nv)', knnsearch(P12_init*Mesh2.verts,Mesh2.verts), ones(Mesh2.nv,1), Mesh2.nv, Mesh1.nv);

end


function [B] = barycentric_coords_w_degenerate(C,V1,V2,V3)
%Finds the barycentric coordinates of triangles and checks for degenerate
%cases. Uses gptoolbox's function barycentric_coordinates.
%C: query point npx3
%V1: corner vertex 1 nvx3
%V2: corner vertex 2 nvx3
%V3: corner vertex 3 nvx3

%first, get the indices where we have degeneracies
B = barycentric_coordinates(C,V1,V2,V3);
id_all = find(vecnorm(V1-V2,2,2) == 0 & (vecnorm(V2-V3,2,2) == 0));
%here, just project to the first index
if(~isempty(id_all))
    B(id_all,1) = 1;
    B(id_all,2:3) = 0;
end
%now, check where only two are equal
id_two = find(vecnorm(V1-V2,2,2) == 0 & (vecnorm(V2-V3,2,2) ~= 0));
%project to V1 and V3
if(~isempty(id_two))
    B(id_two,[1,3]) = barycentric_coordinates(C(id_two,:),V1(id_two,:),V3(id_two,:));
    B(id_two,2)= 0;
end
%check V2 and V3
id_two = find(vecnorm(V1-V2,2,2) ~= 0 & (vecnorm(V2-V3,2,2) == 0));
%project to V1, V3
if(~isempty(id_two))
    B(id_two,[1,3]) = barycentric_coordinates(C(id_two,:),V1(id_two,:),V3(id_two,:));
    B(id_two,2)= 0;
end
%if V1 and V3 are equal
id_two = find(vecnorm(V1-V2,2,2) ~= 0 & (vecnorm(V1-V3,2,2) == 0));
%project to V2, V3
if(~isempty(id_two))
    B(id_two,[1,2]) = barycentric_coordinates(C(id_two,:),V1(id_two,:),V2(id_two,:));
    B(id_two,3)= 0;
end

end


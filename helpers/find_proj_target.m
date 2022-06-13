function [targ, which_tri, bary_tri] = find_proj_target(T,X,Xt)
%find_proj_target
%Inputs:
%       T: target boundary mesh triangulation
%       X: source vertex locations
%       Xt: target vertex locations (on the boundary)
%       return: Xproj: matrix of coordinates of projection targets

% simplices_tri = generate_tri(T, Xt);
% simplices_tri = permute(gather(simplices_tri),[3 1 2]);
% [targ,which_tri,bary_tri] = approximateClosestPointInSimplex(simplices_tri,X,50);
[sqrD,which_tri,targ] = point_mesh_squared_distance(X,Xt,T);
if nargout == 3
   bary_tri = barycentric_coordinates(targ,Xt(T(which_tri,1),:),...
       Xt(T(which_tri,2),:),Xt(T(which_tri,3),:)); 
end
end

function [E, grad,idx] = compute_proj_energy_reverse(X,Xt,T,weights,gamma)
%Computes the projection energy
%Inputs:
%       X: nvx3: source vertex locations. (points we want to project from)
%       Xt: ntx3: target vertex locations (vertices in mesh we want to proj
%       to)
%       T: surface mesh triangulation in target
%       weights: weights per vertex (nvx1)
%       gamma: factor to scale eneryg by
% Projects points in X to mesh (T1,Xt). weights are with respect to X.
%Returns the gradient but with respect to the appropriate target vertices.
[~,I,C] = point_mesh_squared_distance(X,Xt,T);
% get barycentric cooridnates. This function handles degenrate cases.
B = barycentric_coords_w_degenerate(C,Xt(T(I,1),:),Xt(T(I,2),:),Xt(T(I,3),:));
if(isinf(B))
    error('infinite barycentric coords');
end
%TR = triangulation(T,Xt);
%b = cartesianToBarycentric(TR,I,C) ;
%C are the projected points
P = X-C;
%compute the energy
%E = gamma*trace(weights.*(P*P'));
E= gamma*sum(weights.*P.^2,'all');
%compute the gradient with respect to X
grad2 = -2 * gamma *weights.*P;

if nargout > 1
    %get the proper indices
    T1d = T(I,:);
    %compute the gradient with respect to Xt.
    %get the gradient for each vertex
    grad_rep = repmat(grad2,1,1,3);
    grad_rep(:,:,1) = grad_rep(:,:,1).*[B(:,1), B(:,1), B(:,1)];
    grad_rep(:,:,2) = grad_rep(:,:,2).*[B(:,2), B(:,2), B(:,2)];
    grad_rep(:,:,3) = grad_rep(:,:,3).*[B(:,3), B(:,3), B(:,3)];
    %accumulate
    z1 = accumarray(T1d(:),[grad_rep(:,1,1);grad_rep(:,1,2);grad_rep(:,1,3)]);
    z2 = accumarray(T1d(:),[grad_rep(:,2,1);grad_rep(:,2,2);grad_rep(:,2,3)]);
    z3 = accumarray(T1d(:),[grad_rep(:,3,1);grad_rep(:,3,2);grad_rep(:,3,3)]);
    %get the final gradient
    grad = [z1,z2,z3];
    %size is 1xby(max(T1d(:)). Need to get the proper indices.
    idx = unique(T1d(:));
end
end


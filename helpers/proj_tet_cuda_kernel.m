function [P_12, dists, idx, bary] = proj_tet_cuda_kernel(P, T, X)
%Function to project a set of points to a set of tetrahedra
%Inputs: 
%       P: Npxdim: points to project
%       T: Ntx4: tetrahedra indices
%       X: Nvxdim: vertices of the tetrahedra
[P, v0, E_mat, E_dot_mat, E_dot_mat_inv, is_degenerate] = prep_tet_proj_kernel(P, T, X);
[dists,idx,bary] = tet_proj(P,v0,E_mat,E_dot_mat,E_dot_mat_inv,is_degenerate);
bary = reshape(bary,[3,size(bary,1)]);
bary = [1-sum(bary);bary];
idx = double(idx +1); %starts at 0
P_12 = precise_T_to_P_tet(gather([idx,bary']), T, size(gather(X),1));
end
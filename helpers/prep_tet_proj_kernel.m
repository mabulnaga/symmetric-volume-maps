function [P, v0, E_mat, E_dot_mat, E_dot_mat_inv, is_degenerate] = prep_tet_proj_kernel(P, T, X)
%This function prepares the matrices for the tetrahedron projection kernel
%tet_proj.cu
%Inputs: 
%       P: Npxdim: points to project
%       T: Ntx4: tetrahedra indices
%       X: Nvxdim: vertices of the tetrahedra
% Outputs: required matrices to call tet_proj.cu, all GPU Arrays/
%       P: dimxNp: points to project
%       v0: base vertex for each tetrahedron, dimxNt
%       E_mat: matrix for each tetrahedron, v_i - v0. n_txdimx3
%       E_dot_mat: dot product of each combination of E_mat: n_tx3x3
%       E_dot_mat_inv: inverse of E_dot_mat
%       is_degenerate: n_tx1: flag for degenerate tetrahedra

keps =1e-8;
D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
nt = size(T,1);
dim = size(X,2);
dim_simplex = size(T,2);
np = size(P,1);
nv = size(X,1);
%instantiate the matrices we need
%v0 = zeros(nt,dim,'gpuArray');
%is_degenerate = zeros(1,nt,'gpuArray');
%get the base vertex for each tetrahedron.
v0 = gpuArray(X(T(:,1),:));
%compute the simplex basis
Xt = gpu_generate_tets_varin( T,X);
E_matT = pagefun(@mtimes, Xt, D);
E_mat =pagefun(@transpose, E_matT);
E_dot_mat = pagefun(@mtimes, E_mat,E_matT);
%E_dot_mat2 = pagefun(@mtimes, E_mat_T, E_mat);

%check if any of these are degenerate
dets = det3x3(E_dot_mat);
is_degenerate = gpuArray(double(abs(dets) < keps));
E_dot_mat_inv = pagefun(@inv, E_dot_mat);
%is_degenerate = any(all(~isfinite(E_dot_mat_inv)));

%now, properly set the dimensions to obey column-major
P = gpuArray(P');
v0 = v0';
E_mat = permute(E_mat,[2 1 3]);
E_dot_mat = permute(E_dot_mat, [2 1 3]);
E_dot_mat_inv = permute(E_dot_mat_inv, [2 1 3]);

end
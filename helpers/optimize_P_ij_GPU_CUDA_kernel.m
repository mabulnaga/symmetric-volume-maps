function [P_12, P_21] = optimize_P_ij_GPU_CUDA_kernel(X_12, X_21, X_1, X_2, T1, T2, alpha, c_1, c_2, beta)
%Optimize for P_12 and P_21 using a CUDA tet projection kernel. Kernel
%written by Lingxiao Li in "Interactive all-hex meshing via cuboid
%decomposition"
% Inputs:
%        X_12: n1x3 image of mapped vertices of Mesh_1
%        X_21: n2x3 image of mapped vertices of Mesh_2
%        X_1: n1x3 vertices of Mesh_1
%        X_2: n2x3 vertices of Mesh_2
%        T1: nt1x4 connectivity list of tetrahedra of Mesh_1
%        T2: nt2x4 connectivity list of tetrahedra of Mesh_1
%        alpha: scalar: reversibiltiy energy weight
%        c_1: scalar: volume of Mesh_1
%        c_2: scalar: volume of Mesh_2
%        beta: scalar: weight on auxiliary energy.
[A_12, A_21, B_12, B_21] = extract_proj_variables(X_1, X_2, X_12, X_21, c_1, c_2, alpha, beta);
P = B_12;
P_12= proj_tet_cuda_kernel(P, T2, A_12);
P = B_21;
P_21 = proj_tet_cuda_kernel(P, T1, A_21);

end


function [X_tensor, J] = generate_J_tets(T,X,Oinv, use_GPU)
%Wrapper function to compute the Jacobian matrix given a tetrahedral mesh.
%Inputs:
%       T: Ntx4 tetrahdral list
%       X: Nvx3: list of vertices.
%       Oinv: 3x3xNt: tensor containing the original tetrahedron basis.
X_tensor = generate_tets(T,X, use_GPU);
J = compute_J(X_tensor, Oinv, use_GPU);
end


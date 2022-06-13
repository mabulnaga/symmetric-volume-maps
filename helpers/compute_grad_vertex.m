function [grad_X, grad_tet] = compute_grad_vertex(grad_J, D, Oinv, T)
%Functions applies the chain rule to convert a gradient with respect to the
%Jacobian matrix to respect to the vertex positions. We represent the
%Jacobian J as a 3x3 matrix. J = X_i*D(X_{i-1*D})^{-1}. Here, X_i is the
%tetrahedraon at iteration i, X_i is a 3x4 matrix. D is a constant that
%extracts the three basis spanning the tetrahedron.
%Inputs:
%       grad_J: 3x3xNt matrix of gradient with respect to Jacobian
%       Oinv: 3x3xNt matrix of basis vectors of prev tet (X_{i-1}*D)^(-1)
%       T: tetrahedron connectivity list, Ntx3
%Ouputs: grad_X: 3xNv, matrix of gradient per vertex
Nt = size(grad_J,3);
dim = size(grad_J,1);

if(isa(grad_J,'gpuArray'))
    M = pagefun(@mtimes, D, Oinv);
    M = pagefun(@transpose, M);
    grad_tet = pagefun(@mtimes, grad_J, M);
else
  % grad_tet = zeros([size(D),Nt]);
    grad_J = gather(grad_J);
    Oinv = gather(Oinv);
    E = pagemtimes(D,Oinv);
    grad_tet = pagemtimes(grad_J,'none',E,'transpose');
%     for i = 1 : Nt
%         grad_tet(:,:,i) = (grad_J(:,:,i)*(D*Oinv(:,:,i))');
%         %grad_tet(:,:,i) = (D*Oinv(:,:,i))*grad_J(:,:,i);
%     end
end
T = reshape(T',1,[]);
gradMat = reshape(grad_tet,3,[]);
z1 = accumarray(T',gradMat(1,:)')';
z2 = accumarray(T',gradMat(2,:)')';
z3 = accumarray(T',gradMat(3,:)')';
grad_X = [z1',z2',z3'];
end


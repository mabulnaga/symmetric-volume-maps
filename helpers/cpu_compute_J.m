function [ J ] = cpu_compute_J( X,Oinv)
%cpu_compute_J: computes the Jacobian on the CPU.
%   Assumes X, X0 are CPU tensors. Computes the Jacobian as a 3D-array,
%   where the third dimension is the Jacobian of the tetrahedron.
%   Convention: X is 3x4.
D= [-1 -1 -1;1 0 0;0 1 0;0 0 1];
E = pagemtimes(X,D);
J = pagemtimes(E,Oinv);
% for k = 1 : size(X,3)
%    E(:,:,k) =  X(:,:,k) * D;
%    J(:,:,k) = E(:,:,k) * Oinv(:,:,k);
% end

end


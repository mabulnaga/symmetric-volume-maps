function [detX, grad] = determinant_per_vertex(J, T, V)
%DETERMINANT_PER_VERTEX computes the determinant per vertex, and returns
%the gradient
% Inputs: J: 3x3xNt matrix to take the Jacobian
%         T: Ntx4: tetrahedron connectivity list
%         V: Nx3: mesh vertices
%         ID: vertices to take gradient over
% Returns:
%        detX: 1xNt: determinant per tetrahedron
%        grad: Nx3: gradient matrix
N = size(V,1);
ID = 1:N;
D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
%compute the Jacobian matrices
[detX,grad] = det3x3(J);
% now pull to per-vertex
%[detX] = zeros(size(V,1),1);
%V1 = vertexAttachments(triangulation(T,V));
%for i = 1 : length(ID); detX(ID(i)) = min(detJ(V1{ID(i)})); end
% grad per vertex
% only store the ones that we have the ID for

end
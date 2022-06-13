function [alpha] = compute_barycentric_coord(T,X,p)
%Function to compute Barycentric coordinates.
%This ONLY works properly for a n+1 simplex in R^n. For example,
%tetrahedron in R3, triangle in R2. Otherwise, it only enforces the
%constraint sum alpha =1. You can get negative Barycentric coordinates,
%even if inside the tetrahedron.
if(size(T,2) == 3)
    V = [X(:,T(1)), X(:,T(2))];
    C = V - X(:,T(3));
else %if tetrahedron
    V = [X(:,T(1)), X(:,T(2)), X(:,T(3))];
    C = V - X(:,T(4));
end
alpha = linsolve(C,p);
alpha = [alpha; 1-sum(alpha)];
end


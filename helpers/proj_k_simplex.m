function [p_closest,alpha] = proj_k_simplex(V, p)
%Projects to the simplex spanned by columns of V, where V in R^NxK, where
%K are the number of vertices and N is the dimensionality. p is a column
%vector in R^n. We want to find the coordinates \alpha that minimizes 
% \|p-V\alpha\|_2^2 s.t. \sum \alpha = 1.
% %TODO: Add inequlity constraint, \alpha_i >= 0.
if(size(V,2) > size(V,1)+1) %can't project RN to RN+2 simplex
   alpha = zeros(size(p));
   p_closest = p;
   return
end
%make sure p is a column vector. TODO: check if P is a point, or collection
%of points.
if(size(p,2) > 1)
    p = p';
end
one_col = ones(size(V,2),1);
A = [V'*V, -1*one_col; one_col', 0];
B = [V'*p;1];
X = linsolve(A,B);
alpha = X(1:end-1);
p_closest = V*alpha;
end


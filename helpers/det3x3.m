function [det33, grad] =det3x3(A)
%Computes the determinant of a 3x3 matrix
a11 = A(1,1,:);
a12 = A(1,2,:);
a13 = A(1,3,:);
a21 = A(2,1,:);
a22 = A(2,2,:);
a23 = A(2,3,:);
a31 = A(3,1,:);
a32 = A(3,2,:);
a33 = A(3,3,:);

det33 = a11.*a22.*a33 + a12.*a23.*a31 + a13.*a21.*a32 -...
    a13.*a22.*a31 - a12.*a21.*a33 - a11.*a23.*a32;
if(nargout > 1)
    grad = zeros(size(A));
    grad(1,1,:) = a22.*a33-a23.*a32;
    grad(1,2,:) = -(a12.*a33-a13.*a33);
    grad(1,3,:) = a12.*a23-a22.*a13;
    grad(2,1,:) = -(a21.*a33-a31.*a23);
    grad(2,2,:) = a11.*a33-a13.*a31;
    grad(2,3,:) = -(a11.*a23 - a21.*a13);
    grad(3,1,:) = (a21.*a32-a31.*a22);
    grad(3,2,:) = -(a11.*a32-a31.*a12);
    grad(3,3,:) = -(a11.*a22-a12.*a21);
    grad = permute(grad,[2 1 3]);
end
end


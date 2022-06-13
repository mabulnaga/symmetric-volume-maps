function [ T ] = preprocess_flip_volume( T, X )
%PREPROCESS_FLIPVOLUME corrects orientation of tetrahedron volume.
%Inputs:
%       T: Ntx4: tetrahedron connectivity list
%       X: Nx3: N vertices in R^3.
%Outputs:
%       T: Ntx4 tetrahedron with corrected orientation.
%Input a tet, output a tet with all positive signed volumes.

Xvol = volume(X,T);
indices = find(Xvol<0);
for i = 1: length(indices)
    l = T(indices(i),2);
    l2 = T(indices(i),4);
    T(indices(i),2) = l2;
    T(indices(i),4) = l;
end

end


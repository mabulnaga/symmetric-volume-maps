function [O, Oinv, OT, OinvT, OTO] = helper_tensors_SDE(V, use_GPU)
%Computes a few tensors needed to compute the Jacobian matrix and the
%Symmetric Dirichlet Energy. Only works on GPU, and assumes page/tensor
%structure.
%Inputs: matrix of meshes vertices. Size Nx3.
%TODO: Make output variable names more desciptive!
D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
if(use_GPU)
    O = pagefun(@mtimes, V,D);
    Oinv = pagefun(@inv, O);
    OT = pagefun(@transpose, O);
    OinvT = pagefun(@transpose, Oinv);
    OTO = pagefun(@mtimes, OT,O);
else
    O = pagemtimes(V,D);
    OT = pagetranspose(O);
    OTO = pagemtimes(OT,O);
    Oinv = zeros(size(O));
    OinvT = zeros(size(OT));
    for i = 1 : size(O,3)
        Oinv(:,:,i) = inv(O(:,:,i));
        OinvT(:,:,i) = inv(OT(:,:,i));
    end  
end
end


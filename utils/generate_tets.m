function [ tets ] = generate_tets(T,X, use_GPU)
% Function creates a 3D tensor with tetrahedra, size
% tets is size Nx4xNT
% T: size NTx4
% X: size NvxN
% N is the dimension, ex 3 for R^3, or 6 for R^6.
p=1:size(T,1);
p=p';
N = size(X,2);
NT = size(T,1);
if(use_GPU)
    tets = gpuArray(zeros(NT*N,4));
else
    tets = zeros(NT*N,4);
end
tets(:,1) = reshape((X(T(p,1),:))',[],1);
tets(:,2) = reshape((X(T(p,2),:))',[],1);
tets(:,3) = reshape((X(T(p,3),:))',[],1);
tets(:,4) = reshape((X(T(p,4),:))',[],1);
tets = reshape(tets',4,N,[]);
if(use_GPU)
    tets = pagefun(@transpose, tets);
else
    tets = permute(tets,[2 1 3]);
end
end


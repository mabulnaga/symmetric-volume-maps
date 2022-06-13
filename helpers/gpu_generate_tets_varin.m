function [ tets ] = gpu_generate_tets_varin( T,X)
%Generates a tensor where each page is the tetrahedra, represented as a
%R^dimx4 matrix.
% T: size NTx4
% X: size NVxdim
%Output: dimx4xNT.
p=1:size(T,1);
p=p';
tets = gpuArray(zeros(size(T,1)*size(X,2),size(T,2)));
tets(:,1) = reshape((X(T(p,1),:))',[],1); 
tets(:,2) = reshape((X(T(p,2),:))',[],1); 
tets(:,3) = reshape((X(T(p,3),:))',[],1); 
tets(:,4) = reshape((X(T(p,4),:))',[],1); 
tets = reshape(tets',4,size(X,2),[]);
tets = pagefun(@transpose, tets);
end


function [ tets ] = gpu_generate_tets( T,X)
%Generates a tensor with all of the tetrahedra.

p=1:size(T,1);
p=p';
tets = gpuArray(zeros(size(T,1)*3,size(T,2)));
tets(:,1) = reshape((X(T(p,1),:))',[],1); 
tets(:,2) = reshape((X(T(p,2),:))',[],1); 
tets(:,3) = reshape((X(T(p,3),:))',[],1); 
tets(:,4) = reshape((X(T(p,4),:))',[],1); 
tets = reshape(tets',4,3,[]);
tets = pagefun(@transpose, tets);
end


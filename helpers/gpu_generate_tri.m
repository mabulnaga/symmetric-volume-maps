function [ tris ] = gpu_generate_tri( T,X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% tets = zeros(3,4,size(T,1));
% tets = arrayfun(@(x) create_tet(T(x,:),X), 1:size(T,1),'UniformOutput', false);
% tets = reshape(cell2mat(tets),[3,4,size(T,1)]);
% tets = gpuArray(tets);
% tets2 = tets;
p=1:size(T,1);
p=p';
tri = gpuArray(zeros(size(T,1)*size(X,2),size(T,2)));
tri(:,1) = reshape((X(T(p,1),:))',[],1); 
tri(:,2) = reshape((X(T(p,2),:))',[],1); 
tri(:,3) = reshape((X(T(p,3),:))',[],1); 
tris = reshape(tri',3,size(X,2),[]);
tris = pagefun(@transpose, tris);
end


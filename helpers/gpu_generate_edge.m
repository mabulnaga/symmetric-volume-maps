function [ edges ] = gpu_generate_edge( E,X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% tets = zeros(3,4,size(T,1));
% tets = arrayfun(@(x) create_tet(T(x,:),X), 1:size(T,1),'UniformOutput', false);
% tets = reshape(cell2mat(tets),[3,4,size(T,1)]);
% tets = gpuArray(tets);
% tets2 = tets;
% Edge: 2XNe
p=1:size(E,1);
p=p';
edges = gpuArray(zeros(size(E,1)*size(X,2),size(E,2)));
edges(:,1) = reshape((X(E(p,1),:))',[],1); 
edges(:,2) = reshape((X(E(p,2),:))',[],1); 
edges= reshape(edges',2,size(X,2),[]);
edges = pagefun(@transpose, edges);
end


function [P_12] = optimize_P_ij_approx(P_12, X_12, X_21, X_1, X_2, T2, T1s, T2s, alpha, c_1, c_2, beta)
%OPTIMIZE__PIJ_APPROX Uses an approximation to minize P_12 on the CPU.
%Inputs:
% P_12: n1xn2 map from Mesh_1 to Mesh_2
% X_12: n1x3 image of Mesh_1 under the map
% X_21: n2x3 image of Mesh_2 under the map
% X_1: n1x3 Mesh_1 vertices
% X_2: n2x3 Mesh_2 vertices
% T2: ntx4: tetrahdron connectivity of Mesh_2
% T1s: nfx3: boundary face connectivity of Mesh_1
% T2s: nfx3: boundary face connectivity of Mesh_2
% alpha: scalar: reversibility energy weight parameter
% c_1: scalar: total volume of Mesh_1
% c_2: scalar: total volume of Mesh_2
% beta: scalar: parameter on auxiliary energy
% Output:
% P_12: n1xn2: optimized map

% number of barycentric coordinate samples. Increase to improve accuracy
n_sample = 50;
a_1 = (1-alpha)/c_1^2;
b = beta/c_1/c_2;

A_12 = [sqrt(a_1)*X_21, sqrt(b)*X_2];

B_12 = [sqrt(a_1)*X_1, sqrt(b)*X_12];

%minimize for P_12 : projection for each row. P_12: map from M1 to M2,
% constrained to lie on faces of M2. So, this becomes a tet projection
% problem, where we want to project each row of B_12 (vector in R^6) onto each
% tetrahedron of M2, with vertices given by A_12
%loop for each point, over each tetrahedron.
%first, find which vertices are on the boundary of mesh 1.

boundaryVerts1 = unique(T1s(:));
P_12(:) = 0;

batch_size = size(B_12,1);
nbatch =  floor(size(B_12,1)/batch_size);
rem_batch = rem(size(B_12,1),batch_size);
%need support for triangle mesh here too.
bary = zeros(size(B_12,1),4);
which_all = zeros(1,size(B_12,1));
%generate tet simplices once
simplices = generate_tets(T2,A_12,false);
simplices = permute(gather(simplices),[3 1 2]);
for i = 1 : nbatch+1
    if(i <= nbatch)
        inds = (i-1)*batch_size+1:batch_size*i;
    else
        inds = (i-1)*batch_size+1:(i-1)*batch_size+rem_batch;
    end
    query = B_12(inds,:);
    [~,which,bary(inds,:)] = approximateClosestPointInSimplex(simplices,query,n_sample);
    which_all(inds) = which;
end
%need to batch boundary checks
if(~isempty(T2s))
    simplices_tri = generate_tri(T2s,A_12,false);
    simplices_tri = permute(gather(simplices_tri),[3 1 2]);
    [~,which_tri,bary_tri] = closestSimplex(simplices_tri,query);
end

for i = 1 : size(P_12,1)
    if(ismember(i,boundaryVerts1))
        P_12(i,T2s(which_tri(i),:)) = bary_tri(i,:);
    else
        P_12(i,T2(which_all(i),:)) = bary(i,:);
    end
end


end


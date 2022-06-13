function [Mesh1, Mesh2] = prep_mesh_exp(Mesh1,Mesh2)
%Prepares the meshes for experiments. Computes voronoi volumes, areas, and
%tensors for computing the Jacobian.
% Returns:
%         Mesh1, Mesh2: structs contanining mesh parameters.

%[L_1, M1 , star1] = GeometricPrimalLM(Mesh1);
%[L_2, M2, star2] = GeometricPrimalLM(Mesh2);
%barycentric volumes
use_GPU.gpuAvail=false;
V_1 = voronoi_volume(Mesh1.tets, Mesh1.verts);
V_1 = sparse(1:length(V_1), 1:length(V_1), V_1);
V_2 = voronoi_volume(Mesh2.tets, Mesh2.verts);
V_2 = sparse(1:length(V_2), 1:length(V_2), V_2);
% barycentric surface areas
A_1 = voronoi_area(Mesh1.boundary_faces, unique(Mesh1.boundary_faces(:)), Mesh1.verts);
A_2 = voronoi_area(Mesh2.boundary_faces, unique(Mesh2.boundary_faces(:)), Mesh2.verts);
%per tetrahedron volumes
Mesh1.V_1 = V_1;
Mesh2.V_1 = V_2;
Mesh1.A_1 = A_1;
Mesh2.A_1 = A_2;
Mesh1.A_1_normalized = A_1./sum(A_1);
Mesh2.A_1_normalized = A_2./sum(A_2);
%total area
Mesh1.s = sum(Mesh1.A_1);
Mesh2.s = sum(Mesh2.A_1);
%boundary indices
Mesh1.bd_indices = false(length(Mesh1.verts),1);
Mesh1.bd_indices(unique(Mesh1.boundary_faces(:))) = 1;
Mesh2.bd_indices = false(length(Mesh2.verts),1);
Mesh2.bd_indices(unique(Mesh2.boundary_faces(:))) = 1;

% pre-compute tensors for deformation
Mesh1.X0_tensor = generate_tets(Mesh1.tets, Mesh1.verts, false);
% normalized volume
Mesh1.vol_normalized = arrayfun(@(x) abs(tet_volume_signed(Mesh1.X0_tensor(:,:,x))), 1:size(Mesh1.X0_tensor,3));
Mesh1.vol_normalized = Mesh1.vol_normalized/sum(Mesh1.vol_normalized);
Mesh1.A_1_normalized = Mesh1.A_1./Mesh1.s;
% tensors to compute the Jacobian
[Mesh1.O, Mesh1.Oinv, Mesh1.OT, Mesh1.OinvT, Mesh1.OTO] = helper_tensors_SDE(Mesh1.X0_tensor,use_GPU.gpuAvail);
% repeat for Mesh 2
Mesh2.X0_tensor = generate_tets(Mesh2.tets, Mesh2.verts, false);
Mesh2.vol_normalized = arrayfun(@(x) abs(tet_volume_signed(Mesh2.X0_tensor(:,:,x))), 1:size(Mesh2.X0_tensor,3));
Mesh2.vol_normalized = Mesh2.vol_normalized/sum(Mesh2.vol_normalized);
Mesh2.A_1_normalized = Mesh2.A_1/Mesh2.s;
[Mesh2.O, Mesh2.Oinv, Mesh2.OT, Mesh2.OinvT, Mesh2.OTO] = helper_tensors_SDE(Mesh2.X0_tensor, use_GPU.gpuAvail);

end


function [Mesh,landmarks] = prepare_mesh_landmarks(data_dir, generate_landmarks)
% Prepares mesh files for mapping and provides an interactive prompt to
% mark landmarks. Outputs the mesh data and landmarks. Also saves to the
% input data directory.
% If wanting to manually mark landmarks, click iteratively between mesh 1
% and mesh 2. Click 'q' when finished.
%Inputs:
%       data_subdir: full path to folder where meshes are stored. Folder
%       must contain two meshes with .VTK extension
%       generate_landmarks: 0 to skip, 1 for interactive prompt, 2 for closest
%       point matching using k-nearest neighbors.
%Outputs:
%        Mesh: 1x2 cell containing mesh data for mapping
%        landmarks: 1x2 struct containing Euclidean coordinates of marked
%           landmark points.

file_extension = '.vtk';
% read files in the data folder
data_files = extract_dir_files(data_dir,1);
data_files = data_files(arrayfun(@(x) ~contains(lower(x.name),file_extension)~=1,data_files));
% check if no .VTK files found
if(isempty(data_files))
    error('no .VTK meshes found');
end
% sort the meshes: requires natsortfiles (https://www.mathworks.com/matlabcentral/fileexchange/47434)
% [~,ind] = natsortfiles({data_files.name});
ind = 1 : length(data_files);

% make sure there are only two meshes in the data directory. 
if(length(ind)~=2)
    error('incorrect number of meshes found. Requires 2.');
end
% extract these files, create the meshes. Assumes these are all ordered.
% create the meshes
for ii = 1:length(data_files)
    [node, elem] = read_vtk([data_files(ii).folder,'/',data_files(ind(ii)).name]);
    node = node';
    elem = elem';
    mesh_name = data_files(ind(ii)).name;
    Mesh{ii} = prepare_mesh_shape(elem,node,0);
    Mesh{ii}.name = mesh_name;
end
Mesh1 = Mesh{1};
Mesh2 = Mesh{2};
% Generate the landmarks
if(generate_landmarks>0)
    % call the interactive prompt
    if(generate_landmarks == 1)
        [landmarks.p1, landmarks.p2] = create_landmarks(Mesh1,Mesh2);
        % check if we have the correct lengths
        if(size(landmarks.p1,1)~=size(landmarks.p2,1))
            error('inconcistent number of landmarks');
        end
    else
        %knn to automatically pick the closest points
        bd_verts1 = unique(Mesh1.boundary_faces(:));
        bd_verts2 = unique(Mesh2.boundary_faces(:));
        int_verts1 = setdiff(1:length(Mesh1.verts),bd_verts1);
        int_verts2 = setdiff(1:length(Mesh2.verts),bd_verts2);
        if(length(Mesh1.verts) <=Mesh2.verts)
            idx = knnsearch(Mesh1.verts(bd_verts1,:),Mesh2.verts(bd_verts2,:));
            idx_int = knnsearch(Mesh1.verts(int_verts1,:),Mesh2.verts(int_verts2,:));
            %note that some vertices may not have a matching pair.
            landmarks.p2 = Mesh2.verts(bd_verts2,:);
            landmarks.p1 = Mesh1.verts(idx,:);
            landmarks.p2 = [landmarks.p2; Mesh2.verts(idx_int,:)];
            landmarks.p1 = [landmarks.p1; Mesh1.verts(int_verts2,:)];
        else
            idx = knnsearch(Mesh2.verts(bd_verts2,:),Mesh1.verts(bd_verts1,:));
            idx_int = knnsearch(Mesh2.verts(int_verts2,:),Mesh1.verts(int_verts1,:));
            %map the boundary indices appropriately
            %             m1_map = map_bd_indices_to_all(Ts1,Mesh1.boundary_faces);
            %note that some vertices may not have a matching pair.
            landmarks.p1 = Mesh1.verts(bd_verts1,:);
            landmarks.p2 = Mesh2.verts(idx,:);
            landmarks.p1 = [landmarks.p1; Mesh1.verts(int_verts1,:)];
            landmarks.p2 = [landmarks.p2; Mesh2.verts(idx_int,:)];
            %landmarks.p1 = Mesh1.verts; landmarks.p2 = Mesh2.verts;
        end
    end

    %plot and see the landmarks
    f3 = figure;
    f4 = figure;
    figure(f3);
    trisurf(triangulation(Mesh1.boundary_faces,Mesh1.verts));
    hold on
    axis equal;
    title('mesh 1')
    C = hsv(size(landmarks.p1,1));
    scatter3(landmarks.p1(:,1),landmarks.p1(:,2),landmarks.p1(:,3),512,C,'filled');
    figure(f4);
    trisurf(triangulation(Mesh2.boundary_faces,Mesh2.verts));
    title('mesh 2');
    hold on
    axis equal;
    scatter3(landmarks.p2(:,1),landmarks.p2(:,2),landmarks.p2(:,3),512,C,'filled');
    % save the landmarks
    save([data_files(1).folder,'/landmarks.mat'],'landmarks')
else
    landmarks = [];
end
% save the generated meshes
save([data_files(1).folder,'/meshes.mat'],'Mesh')
end

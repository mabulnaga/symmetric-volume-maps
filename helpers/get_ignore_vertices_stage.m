function [ignore_verts, repair_verts] = get_ignore_vertices_stage(stage,Mesh1,Mesh2)
%GET_IGNORE_VERTICES_STAGE Gets which vertices to ignore in the
%optimization depending on the input stage.
%Inputs: stage: 0 (optimization with fixed boundary), 1 (tet inversion), 2
%               optimize all variables
%        Mesh1, Mesh2: structs containing the meshes
%Outputs:
%        ignore_verts: vertices to ignore in the optimization
%        repair_verts: vertices to use in the tet repair stage.

if(stage ==0 || stage == 2)
    if(stage==0)
        ignore_verts{1} = unique(Mesh1.boundary_faces(:));
        ignore_verts{2} = unique(Mesh2.boundary_faces(:));
    else
        ignore_verts{1} = [];
        ignore_verts{2} = [];
    end
    repair_verts{1} = [];
    repair_verts{2} = [];
else
    repair_verts{1} = [];
    repair_verts{2} = [];
    ignore_verts{1} = [];
    ignore_verts{2} = [];
end
end


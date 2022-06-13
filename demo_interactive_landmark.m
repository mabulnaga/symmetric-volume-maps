%% This demo demonstrates using our interactive tool to manually mark landmarks
% and generate the necessary matlab mesh structure for mapping. Assumes
% data is given as .VTK format.
% Instructions: To run the interactive tool, simply alternate clicking
% landmarks on the boundary of Mesh 1 and Mesh 2. Once you are satisfied
% with the selected landmarks, click 'q' to quit. The tool will then save
% the generated landmarks and .mat file.

%% Data
data_path = './data/centaur_0_to_centaur_1_only_VTK/';

%% Map parameters
alpha = 0.5;
gamma = 25;
beta = 5;
tet_uninv_nring = 1;
energy = @ARAP_energy;
lock_bd = false;
%% Mesh preparation

manually_mark_landmarks = 1;
% call function to convert .VTK to mesh data structure. 
[Mesh,landmarks] = prepare_mesh_landmarks(data_path, manually_mark_landmarks);
%% Map!
[X_12, X_21, P_12, P_21, E] = symmetric_volume_map(Mesh, alpha, gamma, beta, landmarks, energy, lock_bd, tet_uninv_nring);

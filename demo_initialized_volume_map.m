%% This demo demonstrates mapping using provided data with an initial volume map,
% % landmarks, and mesh matlab files given. The volume map was initialized
% by first mapping the boundary meshes using the Reversible Harmonic Maps
% method by Danielle Ezuz et al. (2019), then mapping interior vertices to
% the boundary.
clear all; clc; close all;
data_path = './data/centaur_0_to_centaur_1/';
%% Map parameters
alpha = 0.5;
gamma = 25;
beta = 5;
tet_uninv_nring = 1;
energy = @ARAP_energy;
% set this to true to lock the boundary in the first phase of the
% optimization.
lock_bd = true;
%% Load inputs
% load the mesh
Mesh = load([data_path,'meshes.mat']);
Mesh=Mesh.Mesh;
% load the landmarks. Note, if no landmarks are provided, these can be set to
% be an arbitrary vertex:
% landmarks.p1 = Mesh{1}.verts(1,:); landmarks.p2 = Mesh{2}.verts(1,:);
landmarks = load([data_path,'landmarks.mat']);
landmarks = landmarks.landmarks;
% load the initial maps
P12_init= load([data_path,'P12_init.mat']).P12_init;
P21_init= load([data_path,'P21_init.mat']).P21_init;

%% Map!
[X_12, X_21, P_12, P_21, E] = symmetric_volume_map(Mesh, alpha, gamma, beta, landmarks, energy, lock_bd,...
    tet_uninv_nring, P12_init, P21_init);
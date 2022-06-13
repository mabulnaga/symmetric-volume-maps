function [X_12, X_21, P_12, P_21, E_vec] = symmetric_volume_map(Mesh, alpha, gamma, beta, landmarks, energy, lock_bd, tet_uninv_nring, varargin)
%% Symmetric Volume Map algorithm
%Inputs:
% Mesh: 1x2 cell containing input meshes
% alpha: 1x1 reversibiltiy energy parameter
% gamma: 1x1 projection energy parameter
% beta: 1x1 auxiliary energy parameter
% landmarks: 1x2 struct containing landmark vertex coordinates
% energy: function handle to distortion energy
% lock_bd: flag to lock boundary in the first phase of optimization
% tet_uninv_nring: 1x1 nring neighborhood of vertices to use in tet repair
% varargin: pre-initialized maps
%         P_12_init: n1xn2 map from Mesh_1 to Mesh_2, where n1 (n2) are the
%           number of vertices in Mesh 1 (mesh 2).
%         P21_init: n2xn1 initialized map from Mesh_2 to Mesh_1
%Outputs:
% X_12: n1x3 mapped vertices of Mesh_1 (used in paper)
% X_21: n2x3 mapped vertices of Mesh_2 (used in paper)
% P_12: n1xn2 constrained map from Mesh_1 to Mesh_2
% P_21: n2xn1 constrained map from Mesh_2 to Mesh_1
% E_vec: vector of energy values during the optimization.
%% Optimization settings and initial maps
fixBoundary = lock_bd;
%Check if the map is already initialized
if(~isempty(varargin))
    if(length(varargin)~=2)
        error('Incorrect number of pre-initialized maps. Requires 2.')
    end
    P12_init = varargin{1};
    P21_init = varargin{2};
    initialized_map = true;
else
    initialized_map = false;
end
num_times_uninvert_tets = 1;
%% convergence parameters
numIterations = 50;
if(initialized_map)
    num_iterations_min = 10;
else
    num_iterations_min = 40;
end
grad_threshold = 1e-5; 
energy_threshold = 1e-6; 
beta_max = beta;
gamma_max = gamma;
%iteration limits on fmincon.
max_lbfgs_iterations = 1000;
main_lbfgs_iterations = 400;
fmincon_threshold = grad_threshold*2;
%set the optimization distortion energies
func_energy = energy;
energy_tet_repair = func_energy;
E_vec = [];
%%  Initialize the map and pre-compute volumes, areas, and tensors
Mesh1 = Mesh{1};
Mesh2 = Mesh{2};

landmarks_1 = landmarks.p1;
landmarks_2 = landmarks.p2;
% Map initialization
landmarks1 = knnsearch(Mesh1.verts, landmarks_1);
landmarks2 = knnsearch(Mesh2.verts, landmarks_2);
Mesh1.landmarks = landmarks1;
Mesh2.landmarks = landmarks2;
if(~initialized_map)
    [P12_init, P21_init] = initialize_map_landmark(Mesh1, Mesh2, landmarks1, landmarks2);
end
%normalize the vertices so that the meshes have volume 1.
Mesh1 = normalize_mesh_vertices(Mesh1);
Mesh2 = normalize_mesh_vertices(Mesh2);

%Precompute volumes and areas
[Mesh1, Mesh2] = prep_mesh_exp(Mesh1,Mesh2);
c_1 = trace(Mesh1.V_1);
c_2 = trace(Mesh2.V_1);

% initialize optimization variables. Note this is done after normalization
X_1 = Mesh1.verts;
X_2 = Mesh2.verts;
X_12 = P12_init * X_2;
X_21 = P21_init * X_1;
P_12 = P12_init;
P_21 = P21_init;
%% Start mapping!
ct=1;
tStart = tic;
num_times_repair = 0;
% maximum number of stages in the optimization.
max_optimization_count = max(3,(1 + 2*num_times_uninvert_tets)) +1;
count_start = 1;
for optimization_counter = count_start : max_optimization_count
    e_diff = Inf;
    grad_norm = Inf;
    iter = 0;
    E_prev = Inf;
    if(optimization_counter > 1)
        beta = beta_max;
        gamma = gamma_max;
    end
    % get the optimization stage.
    [stage] = get_optimization_stage(fixBoundary, optimization_counter, max_optimization_count, num_times_repair, num_times_uninvert_tets);
    % get the vertices to ignore depending on the stage.
    [ignore_verts, repair_verts] = get_ignore_vertices_stage(stage,Mesh1,Mesh2);
    % check convergence for the current stage.
    while (iter <= numIterations && grad_norm > grad_threshold && e_diff > energy_threshold || iter<num_iterations_min ) %
        iter = iter + 1;
        %increase beta, gamma per iteration for first phase
        if(optimization_counter == 1)
            if(beta < beta_max)
                beta = beta_max/num_iterations_min*iter;
            else
                beta = beta_max;
            end
            if(gamma < gamma_max)
                gamma = gamma_max/num_iterations_min*iter;
            else
                gamma = gamma_max;
            end
        end
        % TET REPAIR STAGE
        if(stage == 1 || stage == 3)
            % break the loop. after doing the uninversion
            iter=numIterations+1;
            %keep a count of how many times we've inverted.
            num_times_repair = num_times_repair + 1;
            nring = tet_uninv_nring;
            if(stage == 1)
                options = optimoptions('fmincon','Display','final','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',max_lbfgs_iterations,...
                    'HessianApproximation','lbfgs','StepTolerance',1e-6,'OptimalityTolerance',fmincon_threshold/10); %prev: both 1e-6. Maybe increase as iterations increase...
                optimize_all_flips = false;
            % post-convergence repair
            elseif(stage == 3)
                options = optimoptions('fmincon','Display','final','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',100,...
                    'HessianApproximation','lbfgs','StepTolerance',1e-6,'OptimalityTolerance',fmincon_threshold/10); %prev: both 1e-6. Maybe increase as iterations increase...
                optimize_all_flips = true;
                nring=0;
                energy_tet_repair = @ARAP_energy_negative_6order;
            end
            % tet repair step.
            [X_12, X_21] = inverted_tet_repair(X_1, X_2, X_12, X_21, P_12, P_21, Mesh1, Mesh2, c_1, c_2, nring, energy_tet_repair, options, optimize_all_flips);
        % OPTIMIZATION STAGE
        else
            % compute the energy
            [E,~,ED,ER,EQ,Ep,Ep_rev] = compute_f_all(Mesh1, Mesh2, [X_12(:);X_21(:)], P_12, P_21, Mesh1.V_1, Mesh2.V_1, X_1, X_2, func_energy,...
                c_1, c_2,  alpha, beta, gamma, Mesh1.vol_normalized,Mesh2.vol_normalized, Mesh1.A_1, Mesh2.A_1, ignore_verts, repair_verts);
            E_vec = [E_vec,E];
            %energy difference to test convergence
            e_diff = abs(E-E_prev);
            E_prev = E;
            fprintf('Iteration no. %d, Total Energy: %0.3f, E_D: %0.3f, E_R: %0.3f, E_Q: %0.3f, E_P: %0.3f, E_P_rev: %0.3f \n\r',ct, E, ED, ER, EQ,Ep, Ep_rev );
            % optimize with LBFGS
            if(iter < 5)
                maxIter = max_lbfgs_iterations;
                optim_tol = fmincon_threshold /10;
            else
                maxIter = main_lbfgs_iterations;
                optim_tol = fmincon_threshold;
            end

            %LBFGS fmincon
            options = optimoptions('fmincon','Display','final','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',maxIter,...
                'HessianApproximation','lbfgs','StepTolerance',1e-6,'OptimalityTolerance',optim_tol); %prev: both 1e-6. Maybe increase as iterations increase...
            [X_min] = fmincon(@(X)compute_f_all(Mesh1, Mesh2, X, P_12, P_21, Mesh1.V_1, Mesh2.V_1, X_1, X_2, func_energy,...
                c_1, c_2,  alpha, beta, gamma, Mesh1.vol_normalized,Mesh2.vol_normalized,Mesh1.A_1,Mesh2.A_1, ignore_verts,repair_verts),[X_12(:);X_21(:)],[],[],[],[],-Inf,Inf,[],options);
            % get the variables, recompute the energy
            X_12 = X_min(1:Mesh1.nv*3,:); X_12 = reshape(X_12,[],3);
            X_21 = X_min(Mesh1.nv*3+1:end,:); X_21 = reshape(X_21,[],3);
            [E,grad,ED,ER,EQ, Ep, Ep_rev] = compute_f_all(Mesh1, Mesh2, [X_12(:);X_21(:)], P_12, P_21, Mesh1.V_1, Mesh2.V_1, X_1, X_2, func_energy,...
                c_1, c_2,  alpha, beta, gamma, Mesh1.vol_normalized,Mesh2.vol_normalized,Mesh1.A_1,Mesh2.A_1, ignore_verts, repair_verts);
            grad_norm = norm(grad(:));
            
            E_vec = [E_vec,E];
            E_prev = E;
            fprintf('Iteration no. %d, Total Energy: %0.3f, E_D: %0.3f, E_R: %0.3f, E_Q: %0.3f, E_P: %0.3f, E_P_rev: %0.3f, \n\r',ct, E, ED, ER, EQ,Ep, Ep_rev);

            %Optimize P
            if(gpuDeviceCount>0)
                P_12old =P_12;
                P_21old = P_21;
                tic
                [P_12, P_21] = optimize_P_ij_GPU_CUDA_kernel(X_12, X_21, X_1, X_2, Mesh1.tets, Mesh2.tets, alpha, c_1, c_2, beta);
                fprintf('projection on GPU\n\r');
                if(fixBoundary) %lock the map on the boundary.
                    P_12(ignore_verts{1},:) = P_12old(ignore_verts{1},:);
                    P_21(ignore_verts{2},:) = P_21old(ignore_verts{2},:);
                end
                toc
            else
                tic
                P_12old = P_12;
                P_21old = P_21;
                P_12 = optimize_P_ij_justin_approx(P_12, X_12, X_21, X_1, X_2, Mesh2.tets, [], [], alpha, c_1, c_2, beta);
                if(fixBoundary) %lock the map on the boundary.
                    P_12(ignore_verts{1},:) = P_12old(ignore_verts{1},:);
                    P_21(ignore_verts{2},:) = P_21old(ignore_verts{2},:);
                end
                toc
            end
            ct = ct + 1;
        end
    end
end
timeElapsed = toc(tStart)
end



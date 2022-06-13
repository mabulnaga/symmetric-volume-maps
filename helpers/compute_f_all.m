function [E, grad, ED, ER, EQ, Ep, Ep_rev] = compute_f_all(Mesh1, Mesh2, X, P_12, P_21, V_1, V_2, X_1, X_2,...
    f, c_1, c_2, alpha, beta, gamma,  weights1, weights2, weights1_proj, weights2_proj, ignore_verts, repair_verts )
%Computes all the energies and places in a format that can be optimized
%using fmincon.
%Inputs:
%       Mesh1: struct containing Mesh1
%       Mesh2: struct containing Mesh2
%       X: 3*(n1+n2)x1: concatenation of X_12, X_21
%       P_12: n1xn2 map from Mesh1 to Mesh2
%       P_21: n2xn1 map from Mesh2 to Mesh1
%       V_1: n1xn1 mass matrix for Mesh1
%       V_2: n2xn2 mass matrix for Mesh2
%       X_1: n1x3: vertices of Mesh1
%       X_2: n2x3 vertices of mesh 2
%       f: function handle to distortion energy
%       c_1: total volume of Mesh1
%       c_2: total volume of mesh2
%       alpha: parameter for reversibility energy
%       beta: parameter for auxiliary energy
%       gamma: parameter for projection energy
%       weights1: 1xnt1: scalar with per-tetrahedron weights for Mesh1
%       weights2: 1xnt2: scalar with per-tetrahedron weights for Mesh2
%       weights1_proj: 1xnf1: scalar with per-boundary face weights for
%           Mesh1
%       weights2_proj: 1xnf2: scalar with per-boundary face weights for
%           Mesh2
%       ignore_verts: 1x2 cell of vertex indices to ignore in the
%           optimization for Mesh1, Mesh2
%       repair_verts: 1x2 cell of vertex indices to use in the tetrahedron
%       repair optimization.
%Ouputs:
%       E: scalar: total energy 
%       grad: 3*(n1+n2)x1: vertex gradient
%       ED, ER, EQ, Ep, Ep_rev: Distortion, reversibility, auxiliary,
%           projection, and reverse projection energies.

use_GPU = false;
if nargin == 19
    ignore_verts{1} = [];
    ignore_verts{2} = [];
    repair_verts{1} = [];
    repair_verts{2} = [];
elseif nargin == 20
    repair_verts{1} = [];
    repair_verts{2} = [];
end
D = [-1 -1 -1;1 0 0;0 1 0;0 0 1];
X1 = X(1:Mesh1.nv*3,:);
X2 = X(Mesh1.nv*3+1:end,:);
%compute the Jacobian matrices
X_12 = reshape(X1,[],3);
[~, J1] = generate_J_tets(Mesh1.tets,X_12,Mesh1.Oinv, use_GPU);
X_21 = reshape(X2,[],3);
[~,J2] = generate_J_tets(Mesh2.tets,X_21,Mesh2.Oinv, use_GPU);

%compute the ARAP energies and gradients
[ED1, gradJ1] = compute_f(J1, f, alpha, weights1);
ED1 = sum(ED1);
[ED2, gradJ2] = compute_f(J2, f, alpha, weights2);
ED2 = sum(ED2);
ED = ED1 + ED2;
%compute gradient w.r.t. vertices
grad_1_d = compute_grad_vertex(gradJ1, D, Mesh1.Oinv, Mesh1.tets);
grad_2_d = compute_grad_vertex(gradJ2, D, Mesh2.Oinv, Mesh2.tets);

%compute the reversibility energy and gradients
[ER, grad_1_r, grad_2_r] = compute_reversibility_energy(X_12, X_21, P_12, P_21, V_1, V_2, X_1, X_2, alpha, c_1, c_2);

%compute the auxiliary energy
if(beta~=0)
    EQ = beta/(c_1*c_2)*compute_aux_energy(P_12, double(X_12), X_2, V_1) + beta/(c_1*c_2)*compute_aux_energy(P_21, double(X_21), X_1, V_2);
    [grad_1_q, grad_2_q] = grad_X_quadsplit(double(X_12), double(X_21), P_12, P_21, V_1, V_2, X_1, X_2, beta, c_1, c_2);
else
    EQ=0; grad_1_q=0; grad_2_q = 0;
end
%compute the projection energy
if(gamma~=0)
    [Ep1, grad_1_p] = projection_energy(X_12,Mesh1, Mesh2, gamma,weights1_proj);
    [Ep2, grad_2_p] = projection_energy(X_21,Mesh2, Mesh1, gamma,weights2_proj);
    Ep = Ep1 + Ep2;
    %now, do the reverse projection
    [Ep_rev1, grad_1_p_rev] = projection_energy_reverse(X_12,Mesh1, Mesh2, 10*gamma,weights2_proj);
    [Ep_rev2, grad_2_p_rev] = projection_energy_reverse(X_21,Mesh2, Mesh1, 10*gamma,weights1_proj);
    Ep_rev = Ep_rev1+Ep_rev2;
else
    Ep = 0;
    grad_1_p = 0;
    grad_2_p = 0;
    Ep_rev = 0;
    grad_1_p_rev = 0;
    grad_2_p_rev = 0;
end
% combine everything
grad1 = grad_1_d + grad_1_r + grad_1_q + grad_1_p + grad_1_p_rev;
grad2 = grad_2_d + grad_2_r + grad_2_q + grad_2_p + grad_2_p_rev;
% kill gradients of vertices to keep fixed.
grad1(ignore_verts{1},:) = 0;
grad2(ignore_verts{2},:) = 0;

%for vertices that we want to focus on inverting only, we set these to only
%have the Dirichlet energy boundary.
grad1(repair_verts{1},:) = grad_1_d(repair_verts{1},:);
grad2(repair_verts{2},:) = grad_2_d(repair_verts{2},:);
grad1 = grad1(:);
grad2 = grad2(:);
grad = [grad1;grad2];
E = ED + ER + EQ + Ep + Ep_rev;

end

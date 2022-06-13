function [E] = compute_aux_energy(P_ij, X_ij, X_j, M_i)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
E = trace((X_ij - P_ij*X_j)'*M_i*(X_ij-P_ij*X_j));
end


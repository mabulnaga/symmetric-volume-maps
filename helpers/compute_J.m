function [J] = compute_J(X, Oinv, use_GPU)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(use_GPU)
    J = gpu_compute_J(X,Oinv);
else
    J = cpu_compute_J(X,Oinv);
end


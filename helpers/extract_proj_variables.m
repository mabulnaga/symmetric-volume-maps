function [A_12, A_21, B_12, B_21] = extract_proj_variables(X_1, X_2, X_12, X_21, c_1, c_2, alpha, beta)
%Prepares the variables in R^6 that are used in the tet projection.

a_1 = (1-alpha)/c_1^2;
a_2 = (1-alpha)/c_2^2;
b = beta/(c_1*c_2);

A_12 = [sqrt(a_1)*X_21, sqrt(b)*X_2];
A_21 = [sqrt(a_2)*X_12, sqrt(b)*X_1]; 

B_12 = [sqrt(a_1)*X_1, sqrt(b)*X_12];
B_21 = [sqrt(a_2)*X_2, sqrt(b)*X_21];
end
function [X_12, X_21] = min_X_ij(P_12, P_21, V_1, V_2, L_1, L_2, X_1, X_2, alpha, beta, c_1, c_2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
X_12 = (alpha/c_1*L_1 + (1-alpha)/c_2^2*P_21'*V_2*P_21 + beta/(c_1*c_2)*V_1)\...
    ((1-alpha)/c_2^2*P_21'*V_2 + beta/(c_1*c_2)*V_1*P_12)*X_2;

X_21 = (alpha/c_2*L_2 + (1-alpha)/c_1^2*P_12'*V_1*P_12 + beta/(c_1*c_2)*V_2)\...
    ((1-alpha)/c_1^2*P_12'*V_1 + beta/(c_1*c_2)*V_2*P_21)*X_1;

end


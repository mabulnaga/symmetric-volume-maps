function [grad_X12, grad_X21] = grad_X_quadsplit(X_12, X_21, P_12, P_21, V_1, V_2, X_1, X_2,  beta, c_1, c_2)
%Computes the gradient for X_12 or X_21 in the quadratic splitting term
if(beta ~=0)
    grad_X12 = (beta)/(c_1*c_2)*2*V_1*(X_12-P_12*X_2);
    grad_X21 = (beta)/(c_1*c_2)*2*V_2*(X_21-P_21*X_1);
else
    grad_X12 = 0;
    grad_X21 = 0;
end
end


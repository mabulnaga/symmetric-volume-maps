function [E,grad_X12, grad_X21, ER1, ER2] = compute_reversibility_energy(X_12, X_21, P_12, P_21, V_1, V_2, X_1, X_2, alpha, c_1, c_2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(alpha~=1)
    ER1 = ((1-alpha)/(c_1)^2)*trace((P_12*X_21-X_1)'*V_1*(P_12*X_21-X_1));
    ER2 = ((1-alpha)/(c_2)^2)*trace((P_21*X_12-X_2)'*V_2*(P_21*X_12-X_2));
    E = ER1 + ER2;
    if nargout > 1
        [grad_X12, grad_X21] = grad_X_reversibility(X_12, X_21, P_12, P_21, V_1, V_2, X_1, X_2, alpha,c_1, c_2);
    end
else
    E = 0;
    grad_X12 = 0;
    grad_X21 = 0;
end
end


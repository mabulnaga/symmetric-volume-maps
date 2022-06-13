function [grad_X12, grad_X21] = grad_X_reversibility(X_12, X_21, P_12, P_21, V_1, V_2, X_1, X_2, alpha,c_1, c_2)
%Computes the gradient for X_12 or X_21 in the reversibility term
grad_X12 = (1-alpha)/c_2^2*(2*P_21'*V_2*(P_21*X_12-X_2));
grad_X21 = (1-alpha)/c_1^2*(2*P_12'*V_1*(P_12*X_21-X_1));


end


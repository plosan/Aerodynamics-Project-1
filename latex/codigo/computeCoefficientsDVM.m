function [Cl, Cm] = computeCoefficientsDVM(U_inf, chord, alpha, x_ref, Gamma, vortex)
% Computation of Cl and Cm_ref coefficients given physical data
Cl = 2*sum(Gamma)/(U_inf*chord); % Lift coefficient
Cm = -2*cos(alpha)*(Gamma')*(vortex(:,1)-x_ref)/(U_inf*chord^2); % Moment coefficient
end


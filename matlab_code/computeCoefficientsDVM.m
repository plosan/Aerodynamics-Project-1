function [Cl, Cm] = computeCoefficientsDVM(N, U_inf, chord, alpha, x_ref, Gamma, vortex)


Cl = 2*sum(Gamma)/(U_inf*chord);

Cm = 0;
for i = 1:N
    Cm = Cm + Gamma(i)*(vortex(i) - x_ref)*cos(alpha);
end
Cm = -2*Cm/(U_inf*chord^2);

end


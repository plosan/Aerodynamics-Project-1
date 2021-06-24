function [A0, A1, A2] = computeACoefficients(f, p, alpha)
% Computation of A0, A1 and A2 coefficients for TAT
theta_p = acos(1-2*p);
fun0 = @(theta) (cos(theta)+2*p-1);                 %Function to compute A0
fun1 = @(theta) (cos(theta)+2*p-1).*cos(theta);     %Function to compute A1
fun2 = @(theta) (cos(theta)+2*p-1).*cos(2*theta);   %Function to compute A2
A0 = integral(fun0, theta_p, pi)/(1-p)^2;   
A1 = integral(fun1, theta_p, pi)/(1-p)^2;
A2 = integral(fun2, theta_p, pi)/(1-p)^2;
if p ~= 0   % Prevent division by 0
    A0 = integral(fun0, 0, theta_p)/p^2 + A0;
    A1 = integral(fun1, 0, theta_p)/p^2 + A1;
    A2 = integral(fun2, 0, theta_p)/p^2 + A2;
end
A0 = alpha - (f/pi)*A0; % Coefficient A0
A1 = (2*f/pi)*A1;       % Coefficient A1
A2 = (2*f/pi)*A2;       % Coefficient A2
end
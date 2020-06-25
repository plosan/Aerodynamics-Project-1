function [A0, A1, A2] = computeACoefficients(f, p, alpha)

theta_p = acos(1 - 2*p);
fun0 = @(theta) (cos(theta)+2*p-1);
% A0 = alpha - (f/pi)*(integral(fun0, 0, theta_p)/p^2 + integral(fun0, theta_p, pi)/(1-p)^2);

fun1 = @(theta) (cos(theta)+2*p-1).*cos(theta);
% A1 = (2*f/pi)*(integral(fun1, 0, theta_p)/p^2 + integral(fun1, theta_p, pi)/(1-p)^2);

fun2 = @(theta) (cos(theta)+2*p-1).*cos(2*theta);
% A2 = (2*f/pi)*(integral(fun2, 0, theta_p)/p^2 + integral(fun2, theta_p, pi)/(1-p)^2);

A0 = integral(fun0, theta_p, pi)/(1-p)^2;
A1 = integral(fun1, theta_p, pi)/(1-p)^2;
A2 = integral(fun2, theta_p, pi)/(1-p)^2;
if p ~= 0
    A0 = integral(fun0, 0, theta_p)/p^2 + A0;
    A1 = integral(fun1, 0, theta_p)/p^2 + A1;
    A2 = integral(fun2, 0, theta_p)/p^2 + A2;
end
A0 = alpha - (f/pi)*A0;
A1 = (2*f/pi)*A1;
A2 = (2*f/pi)*A2;


end
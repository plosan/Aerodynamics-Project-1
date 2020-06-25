clear;
close all;
clc;

%% 1. DATA
% Airfoil
NACA = 2408;
x_flap = 0;
eta = 0;
eta = eta*pi/180;

% Numerical data
distribution = "fullcosine";
N = 200;

% Physical data
U_inf = 1;
chord = 1;
x_ref = 1/4;

% Airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)

%% 2. COMPUTATION FOR DIFFERENT ANGLES OF ATTACK

alpha_step = 1e-2;
alpha_lim = 30;
alpha = -alpha_lim:alpha_step:alpha_lim;
alpha = alpha*pi/180;
Cl_DVM = zeros(1, length(alpha));
Cm0_DVM = zeros(1, length(alpha));

[x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N, distribution); %#ok<ASGLU>

% Discrete Vortex Method
for i = 1:length(alpha)
    Gamma = computeCirculation(U_inf, alpha(i), vortex, node, n_vec, N);  
    [Cl_DVM(i), Cm0_DVM(i)] = computeCoefficientsDVM(N, U_inf, chord, alpha(i), x_ref, Gamma, vortex);
end


%% 3. COMPUTE Cl_alpha, alpha_l0, Cm0

% Angle of 0 lift (alpha_l0)
i = 0;
found = 0;
while (i <= length(alpha)) && (found == 0)
    i = i + 1;
    if Cl_DVM(i) >= 0
        found = 1;
    end
end
x = Cl_DVM(i-1)/(Cl_DVM(i-1)-Cl_DVM(i));                % Linear interpolation          
alpha_l0 = x*alpha(i)+(1-x)*alpha(i-1);     % Angle of 0 lift
alpha_l0 = alpha_l0*180/pi;

% Compute Cl_alpha
lin_lim = 10*pi/180;
found = 0;
i = 1;
i1 = 1;
i2 = 1;
Cl_1 = 0;
Cl_2 = 0;
Cm0_1 = 0;
Cm0_2 = 0;
while (i <= length(alpha)) && (found < 2)
    if alpha(i) == -lin_lim
        found = found + 1;
        i1 = i;
        Cl_1 = Cl_DVM(i);
    elseif alpha(i) == lin_lim
        found = found + 1;
        i2 = i;
    end
    i = i + 1;
end
Cl_alpha = (Cl_DVM(i2) - Cl_DVM(i1))/(2*lin_lim*180/pi);
% err_Cl_alpha = 100*abs(Cl_alpha*180/pi - 2*pi)/(2*pi);

% Compute Cm0
Cm0 = sum(Cm0_DVM(i1:1:i2))/length(Cm0_DVM(i1:1:i2));
min_Cm0 = min(Cm0_DVM);
max_Cm0 = max(Cm0_DVM(i1), Cm0_DVM(i2));
max_variation = abs(max_Cm0 - min_Cm0);

% Cl equation
Cl_equation = sprintf("Cl = %.4f*alpha %.4f", Cl_alpha, Cl_alpha*alpha_l0);

% Print coefficients
fprintf("%10s%15s%15s%15s\n", "i", "alpha", "Cl", "Cm0");
fprintf("%10d%15.3f%15.3f%15.3f\n", i1, alpha(i1)*180/pi, Cl_DVM(i1), Cm0_DVM(i1));
fprintf("%10d%15.3f%15.3f%15.3f\n", i2, alpha(i2)*180/pi, Cl_DVM(i2), Cm0_DVM(i2));
fprintf("%15s = %.9f %s\n", "alpha_l0", alpha_l0, "deg");
fprintf("%15s = %.9f %s\n", "Cl_alpha", Cl_alpha, "1/deg");
fprintf("%15s = %.9f %s\n", "Cl_alpha", Cl_alpha*180/pi, "1/rad");
% fprintf("%15s = %.9f %s\n", "err_Cl_alpha", err_Cl_alpha, "%");
fprintf("%15s = %.9f \n", "Cm0", Cm0);
fprintf("%15s = %.9f \n", "max_variation", max_variation);
fprintf("%30s\n", Cl_equation);


%% 4. PLOT OF LIFT COEFFICIENT AND FREE MOMENT COEFFICIENT
plotDVMCoefficients(alpha, Cl_DVM, Cm0_DVM, i1, i2);







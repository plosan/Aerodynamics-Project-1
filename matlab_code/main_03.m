clear;
close all;
clc;

%% 1. DATA
% Airfoil
NACA = 2408;        % Airfoil
x_flap = 0;         % Flap location
eta = 0*pi/180;     % Flap deflection angle

% Physical data
U_inf = 1;      % Free stream velocity
chord = 1;      % Airfoil chord length
x_ref = 1/4;    % Reference point for moment computation (a.c.)

% Airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)

% Numerical data
N = 200;                        % Number of points
distribution = "fullcosine";    % Type of distribution

%% 2. COMPUTATION FOR DIFFERENT ANGLES OF ATTACK WITHOUT FLAP
alpha_step = 1e-2;  % Step for alpha vector              
alpha_lim = 30;     % Limit for alpha vector
alpha = -alpha_lim:alpha_step:alpha_lim;    % Alpha vector
alpha = alpha*pi/180;                       % Conversion to radians
Cl_DVM = zeros(1, length(alpha));   % Computed lift coefficient
Cm0_DVM = zeros(1, length(alpha));  % Computed free moment coefficient
% Geometry computation
[x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N, distribution); %#ok<ASGLU>
% Discrete Vortex Method
for i = 1:length(alpha)
    Gamma = computeCirculation(U_inf, alpha(i), vortex, node, n_vec, N);  
    [Cl_DVM(i), Cm0_DVM(i)] = computeCoefficientsDVM(U_inf, chord, alpha(i), x_ref, Gamma, vortex);
end

%% 3. COMPUTE alpha_l0, Cl_alpha WITHOUT FLAP
% Computation of 0 lift angle (alpha_l0)
i = 0;  % Initial index for search
found = 0;  % Binary variable to end search
while (i <= length(alpha)) && (found == 0)
    i = i + 1;
    if Cl_DVM(i) >= 0
        found = 1;
    end
end
x = Cl_DVM(i-1)/(Cl_DVM(i-1)-Cl_DVM(i)); % Linear interpolation          
alpha_l0 = x*alpha(i)+(1-x)*alpha(i-1);  % Angle of 0 lift
alpha_l0 = alpha_l0*180/pi; % Conversion to degrees

% Compute Cl_alpha
lin_lim = 10*pi/180; % Linear limit
found = 0; % Binary variable to end search
i = 1;  % Initial index for search
i1 = 1; % Index for lower limit of linear range
i2 = 1; % Index for upper limit of linear range
while (i <= length(alpha)) && (found < 2)
    if alpha(i) == -lin_lim
        found = found + 1;
        i1 = i;
    elseif alpha(i) == lin_lim
        found = found + 1;
        i2 = i;
    end
    i = i + 1;
end
Cl_alpha = (Cl_DVM(i2) - Cl_DVM(i1))/(2*lin_lim);

%% 4. COMPUTATION FOR SEVERAL FLAP-TO CHORD RATIOS 
% Angles of attack for which flap efficiency factor will be computed
alpha_flap = [0]; % Angle of attack for airfoil with flap deflection
% alpha_flap is a vector, so it allows computing flap efficiency factor for
% several angles of attack
alpha_flap = alpha_flap*pi/180; 

% Get lift coefficient without flap from the previous computations
Cl_no_flap = zeros(1, length(alpha_flap));
i = 0;  % Search index
current = 1;    % Current index in alpha_flap
while (i <= length(alpha)-1) && (current <= length(alpha_flap))
    i = i + 1;
    if alpha(i) == alpha_flap(current)
        Cl_no_flap(current) = Cl_DVM(i);
        current = current + 1;
    end
end

% Flap deflection angle
eta = 10;
eta = eta*pi/180;

% Flap chord ratios and flap hinge position
E = 0:0.05:0.4;
xh = chord*(1-E);

% Computation of Cl for airfoil with deflected flap
Cl_flap = zeros(length(xh), length(alpha_flap)); % C
for i = 1:length(xh)
    [x, z, vortex, node, c, n_vec, t_vec] = ...
        computeGeometry(f, p, chord, xh(i), eta, N, distribution);
    % Discrete Vortex Method
    for j = 1:length(alpha_flap) % Compute Cl and Cm0 for each alpha
        Gamma = computeCirculation(U_inf, alpha_flap(j), vortex, node, n_vec, N);  
        [Cl_flap(i,j), Cm0_DVM(i,j)] = ...
            computeCoefficientsDVM(U_inf, chord, alpha_flap(j), x_ref, Gamma, vortex);
    end
end

% Computation of flap effectivenes DVM
flap_eff_DVM = zeros(length(xh), length(alpha_flap));
for i = 1:length(xh)
    flap_eff_DVM(i,:) = (Cl_no_flap - Cl_flap(i,:))/(Cl_alpha*eta);
end

% Computation of corrected flap effectivenes
factor = 0.8;
flap_eff_corrected = factor*flap_eff_DVM;

% Experimental flap effectivenes (from 
flap_eff_exp = -[0, 0.18, 0.29, 0.38, 0.46, 0.54, 0.60, 0.65, 0.70];

%% 5. COMPUTE ERROR
err_DVM = abs((flap_eff_DVM - flap_eff_exp)./flap_eff_exp);
err_corrected = abs((err_DVM - flap_eff_exp)./flap_eff_exp);


for i = 1:length(E)
    fprintf("$%5.2f$", E(i));
    fprintf("%2s&%2s", "", "");
    fprintf("$%8.2f$", -flap_eff_DVM(i));
    fprintf("%2s&%2s", "", "");
    fprintf("$%8.2f$", -flap_eff_corrected(i));
    fprintf("%2s&%2s", "", "");
    fprintf("$%8.2f$", -flap_eff_exp(i));
    fprintf("%2s&%2s", "", "");
    err_DVM = 0;
    err_corrected = 0;
    if flap_eff_exp(i) ~= 0
        err_DVM = 100*abs((flap_eff_DVM(i)-flap_eff_exp(i))/flap_eff_exp(i));
        err_corrected = 100*abs((flap_eff_corrected(i) - flap_eff_exp(i))./flap_eff_exp(i));    
    end
    fprintf("$%8.2f$", err_DVM);
    fprintf("%2s&%2s", "", "");
    fprintf("$%8.2f$", err_corrected);
    fprintf("\t\\\\\n");
end


%% 6. PLOT
figure(1);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Factor de eficiencia del flap $\left( \partial \alpha_{l0} / \partial \eta \right)$}");
for j = 1:length(alpha_flap)
    plot(E, -flap_eff_DVM(:,j), 'b');
end
for j = 1:length(alpha_flap)
    plot(E, -flap_eff_corrected(:,j), 'r');
end
plot(E, -flap_eff_exp, 'g');
scatter(E, -flap_eff_exp, 20, 'g', 'filled');
for j = 1:length(alpha_flap)    
    scatter(E, -flap_eff_DVM(:,j), 20, 'b', 'filled');
end
for j = 1:length(alpha_flap)    
    scatter(E, -flap_eff_corrected(:,j), 20, 'r', 'filled');
end
xlabel("\emph{Flap-chord ratio}");
ylabel("Factor de eficiencia del flap");
set(gcf, 'units', 'centimeters', 'position', [18,1,18,12]);
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.2f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
legend("Factor de eficiencia DVM", ...
    sprintf("Factor de eficiencia DVM corregido $( %.2f )$", factor), ...
    "Factor de eficiencia experimental", ...
    "Location", "Northwest");
grid on;
grid minor;
box on;
hold off;








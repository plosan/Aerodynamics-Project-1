clear;
close all;
clc;

%% 1. DATA
% Airfoil properties
x_flap = 0;         % Flap location
eta = 0*pi/180;     % Flap deflection angle

% Numerical data
N = 200;                    % Number of panels
distribution = "uniform";   % Type of discretization

% Physical data
U_inf = 1;      % Free stream velocity
chord = 1;      % Airfoil chord length
x_ref = 0.25;   % Reference point for moment computation (a.c.)

% Airfoil geometric parameters
f = 0:0.01:0.06;    % Max camber
p = 0.1:0.1:0.6;    % Max camber position


%% 2. COMPUTATION FOR DIFFERENT ANGLES OF ATTACK
alpha_step = 1e-1;  % Step for alpha vector
alpha_lim = 10;     % Limit for alpha vector
alpha = -alpha_lim:alpha_step:alpha_lim;    % Alpha vector
alpha = alpha*pi/180;                       % Conversion to radians
alpha_l0 = zeros(length(p), length(f)); % Angle of 0 lift matrix
Cm0 = zeros(length(p), length(f));      % Free moment coefficient matrix
% Discrete Vortex Method
for i = 1:length(p) % Compute for each p
    for j = 1:length(f) % Compute for each f
        % Compute geometry
        [x, z, vortex, node, c, n_vec, t_vec] = ...
            computeGeometry(f(j), p(i), chord, x_flap, eta, N, distribution);
        % Lift coefficient and free moment
        Cl_DVM = zeros(1, length(alpha));
        Cm0_DVM = zeros(1, length(alpha));
        for k = 1:length(alpha) % Compute for each alpha
            Gamma = computeCirculation(U_inf, alpha(k), vortex, node, n_vec, N);
            [Cl_DVM(k), Cm0_DVM(k)] = ...
                computeCoefficientsDVM(U_inf, chord, alpha(k), x_ref, Gamma, vortex);
        end
        % Compute alpha_l0 (radians) and free moment coefficient
        alpha_l0(i,j) = (Cl_DVM(end)*alpha(1) - Cl_DVM(1)*alpha(end))/(Cl_DVM(end) - Cl_DVM(1));
        Cm0(i,j) = sum(Cm0_DVM)/length(Cm0_DVM);
    end
end
alpha_l0 = alpha_l0*180/pi; % Conversion to degrees

%% 3. PLOTS
% Color vector
color = ['b', 'k', 'r', 'g', 'y', 'c', 'm'];

% Legend strings
legend_str = cell(1,length(f));
for j = 1:length(f)
    legend_str(j) = {sprintf("$f = %.2f$", f(j))};
end

% Plot p vs. alpha_l0, for each f
figure(1);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{\'Angulo de sustentaci\'on nula $\left( \alpha_{l0} \right)$}");
for j = 1:length(f)
    plot(p, alpha_l0(:, j), color(j));
end
for j = 1:length(f)
    scatter(p, alpha_l0(:, j), 20, color(j), 'filled');
end
xlabel("Posici\'on de m\'axima combadura");
ylabel("\'Angulo de sustentaci\'on nula $\left( \mathrm{deg} \right)$");
xticks(p);
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.0f'));
set(gcf, 'units', 'centimeters', 'position', [0,1,18,12]);
legend(legend_str, 'Location', 'Southwest');
grid on;
grid minor;
box on;
hold off;

% Plot p vs. Cm0, for each f
figure(2);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Coeficiente de momento libre $\left( C_{m0} \right)$}");
for j = 1:length(f)
    plot(p, Cm0(:, j), color(j));
end
for j = 1:length(f)
    scatter(p, Cm0(:, j), 20, color(j), 'filled');
end
xlabel("Posici\'on de m\'axima combadura");
ylabel("Coeficiente de momento libre");
xticks(p);
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
set(gcf, 'units', 'centimeters', 'position', [18,1,18,12]);
legend(legend_str, 'Location', 'Southwest');
grid on;
grid minor;
box on;
hold off;








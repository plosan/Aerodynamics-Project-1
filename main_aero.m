clear;
close all;
clc;

% Airfoil
NACA = 2408;
x_flap = 0.7;
eta = 10;

eta = eta*pi/180;

% Numerical data
N = 20;
distribution = "uniform";

% Physical data
U_inf = 1;
chord = 1;
alpha = 4;
alpha = alpha*pi/180;
cos_alpha = cos(alpha);
sin_alpha = sin(alpha);

[x, z, vortex, node, c, n_vec, t_vec, Gamma] = ...
    computeCirculation(NACA, U_inf, chord, x_flap, eta, alpha, N, distribution);

% Compute Cl
sum_gamma = 0;
for i = 1:length(Gamma)
    sum_gamma = sum_gamma + Gamma(i);
end

Cl = 2*sum_gamma/(chord*U_inf);


% Compute delta_Cpi
delta_Cp = (2/U_inf)*(Gamma./c);
% Compute C_mLE
x_ref = 0;
sum_cm = 0;
for i = 1:N-1
    sum_cm = sum_cm + Gamma(i)*(vortex(i) - x_ref)*cos(alpha);
end
Cm_ref = -2*sum_cm/(U_inf*chord^2);


fprintf("%10s = %8.5f\n", "Cl", Cl);
fprintf("%10s = %8.5f\n", "Cm_LE", Cm_ref);

% Plot Airfoil
str = sprintf("\\textbf{NACA %d}", NACA);
figure(1);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title(str);
plot(x, z, 'r');
scatter(x, z, 10, 'r', 'filled');
scatter(vortex(:,1), vortex(:,2), 10, 'g', 'filled');
scatter(node(:,1), node(:,2), 10, 'b', 'filled');
quiver(x(2:end), z(2:end), t_vec(:,1)', t_vec(:,2)');
quiver(node(:,1), node(:,2), n_vec(:,1), n_vec(:,2));
axis equal;
xlim([0 1]);
ylim([-0.5 0.5]);
xlabel("$x / c$");
ylabel("$z / c$");
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.1f'));
grid on;
box on;
set(gcf, 'units', 'centimeters', 'position', [25,1,20,20]);
legend("Camber line", "Camber line points", "Vortex nodes", ...
    "Control nodes", "Tangent vectors", "Normal vectors");
hold off;

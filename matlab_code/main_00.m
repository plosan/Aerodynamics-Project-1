clear;
close all;
clc;

%% 1. DATA

% Airfoil
NACA = 0012;
x_flap = 0;
eta = 0;
eta = eta*pi/180;

% Numerical data
N = 5;                      % Number of panels
distribution = "uniform";   % Type of discretization

% Physical data
U_inf = 1;
chord = 1;
alpha = 5;
alpha = alpha*pi/180;
cos_alpha = cos(alpha);
sin_alpha = sin(alpha);
x_ref = 0;

% Airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)
t = mod(NACA, 100)/100;             % Thickness (percent of chord)


%% 2. DISCRETE VORTEX METHOD 
[x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N, distribution);
[Gamma] = computeCirculation(U_inf, alpha, vortex, node, n_vec, N);
[Cl_DVM, Cm_DVM] = computeCoefficientsDVM(N, U_inf, chord, alpha, x_ref, Gamma, vortex);








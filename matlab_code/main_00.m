clear;
close all;
clc;

%% 1. DATA
% Airfoil
NACA = 0012;        % Airfoil
x_flap = 0;         % Flap location
eta = 0*pi/180;     % Flap deflection angle

% Numerical data
N = 5;                      % Number of panels
distribution = "uniform";   % Type of discretization

% Physical data
U_inf = 1;          % Free stream velocity
chord = 1;          % Airfoil chord length
alpha = 5*pi/180;   % Angle of attack
x_ref = 0;          % Reference point for moment computation 

% Airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)
t = mod(NACA, 100)/100;             % Thickness (percent of chord)

%% 2. DISCRETE VORTEX METHOD 
[x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N, distribution);
[Gamma] = computeCirculation(U_inf, alpha, vortex, node, n_vec, N);
[Cl_DVM, Cm_DVM] = computeCoefficientsDVM(U_inf, chord, alpha, x_ref, Gamma, vortex);








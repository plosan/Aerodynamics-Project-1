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

% Physical data
U_inf = 1;
chord = 1;
alpha = 4*pi/180;

% Airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)


%% 2. THIN AIRFOIL THEORY COMPUTATIONS
[A0, A1, A2] = computeACoefficients(f, p, alpha);
Cl_TAT = (2*A0+A1)*pi;
Cm_TAT = -Cl_TAT/4 + (A2-A1)*pi/4;

%% 3. DISCRETE VORTEX METHOD CONVERGENCE
N = 1:1:200;
Cl_DVM = zeros(1, length(N));
Cm_DVM = zeros(1, length(N));

for i = 1:length(N)
    % Discrete Vortex Method
    [x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N(i), distribution);
    [Gamma] = computeCirculation(U_inf, alpha, vortex, node, n_vec, N(i));  
    [Cl_DVM(i), Cm_DVM(i)] = computeCoefficientsDVM(N(i), U_inf, chord, alpha, 0, Gamma, vortex);
end

%% 4. PLOT ERROR
plotError(N, Cl_TAT, Cl_DVM, Cm_TAT, Cm_DVM);










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

N = 10;

[x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N, distribution);

plotAirfoilDVM(NACA, x, z, vortex, node, n_vec, t_vec);
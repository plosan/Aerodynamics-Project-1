%% 1. DATA
% Airfoil
NACA = 2408;        % Airfoil
x_flap = 0;         % Flap location
eta = 0*pi/180;     % Flap deflection angle

% Numerical data
distribution = "fullcosine";

% Physical data
U_inf = 1;          % Free stream velocity
chord = 1;          % Airfoil chord length
alpha = 4*pi/180;   % Angle of attack

% Airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)

%% 2. THIN AIRFOIL THEORY COMPUTATIONS
[A0, A1, A2] = computeACoefficients(f, p, alpha); %A0, A1, A2 coefficients
Cl_TAT = (2*A0+A1)*pi;  % Lift coefficient
Cm_TAT = -Cl_TAT/4 + (A2-A1)*pi/4;  % Leading edge moment coefficient

%% 3. DISCRETE VORTEX METHOD CONVERGENCE
N = 1:1:200;    % Number of panels
Cl_DVM = zeros(1, length(N));   % Computed lift coefficient
Cm_DVM = zeros(1, length(N));   % Computed leading edge moment coefficient

for i = 1:length(N)
    % Discrete Vortex Method
    [x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N(i), distribution);
    [Gamma] = computeCirculation(U_inf, alpha, vortex, node, n_vec, N(i));  
    [Cl_DVM(i), Cm_DVM(i)] = computeCoefficientsDVM(U_inf, chord, alpha, 0, Gamma, vortex);
end










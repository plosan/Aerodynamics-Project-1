%% 1. DATA
% Airfoil
NACA = 2408;        % Airfoil
x_flap = 0;         % Flap location
eta = 0*pi/180;     % Flap deflection angle

% Numerical data
N = 200;                        % Number of panels
distribution = "fullcosine";    % Type of discretization

% Physical data
U_inf = 1;      % Free stream velocity
chord = 1;      % Airfoil chord length   
x_ref = 1/4;    % Reference point for moment computation (a.c.)

% Airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)

%% 2. COMPUTATION FOR DIFFERENT ANGLES OF ATTACK
alpha_step = 1e-2;  % Step for alpha vector
alpha_lim = 30;     % Limit for alpha vector
alpha = -alpha_lim:alpha_step:alpha_lim;    % Alpha vector
alpha = alpha*pi/180;                       % Conversion to radians
Cl_DVM = zeros(1, length(alpha));           % Computed lift coefficient
Cm0_DVM = zeros(1, length(alpha));          % Computed free moment coefficient
% Geometry computation
[x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N, distribution); %#ok<ASGLU>
% Discrete Vortex Method, compute Gamma, Cl and Cm0 for each alpha
for i = 1:length(alpha)
    Gamma = computeCirculation(U_inf, alpha(i), vortex, node, n_vec, N);  
    [Cl_DVM(i), Cm0_DVM(i)] = computeCoefficientsDVM(U_inf, chord, alpha(i), x_ref, Gamma, vortex);
end

%% 3. COMPUTE Cl_alpha, alpha_l0, Cm0
% Computation of 0 lift angle (alpha_l0)
i = 0; % Initial index for search
found = 0; % Binary variable to end search
while (i <= length(alpha)) && (found == 0)
    i = i + 1;
    if Cl_DVM(i) >= 0
        found = 1;
    end
end
x = Cl_DVM(i-1)/(Cl_DVM(i-1)-Cl_DVM(i)); % Linear interpolation          
alpha_l0 = x*alpha(i)+(1-x)*alpha(i-1);  % Angle of 0 lift
alpha_l0 = alpha_l0*180/pi; % Conversion to degrees

% Find indexs for limits of linear range
lin_lim = 10*pi/180; % Limit of Cl linear range
i = 1; % Initial index for search
i1 = 1; % Index for lower limit of linear range
i2 = 1; % Index for upper limit of linear range
found = 0; % Binary variable to end search  
% Index search
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

% Compute lift slope (1/deg)
Cl_alpha = (Cl_DVM(i2) - Cl_DVM(i1))/(2*lin_lim*180/pi);

% Compute Cm0
Cm0 = sum(Cm0_DVM(i1:1:i2))/length(Cm0_DVM(i1:1:i2));
min_Cm0 = min(Cm0_DVM); % Minimum Cm0
max_Cm0 = max(Cm0_DVM(i1), Cm0_DVM(i2)); % Max Cm0 
% Max relative variation of Cm0 wrt minimum Cm0
max_variation = 100*abs((max_Cm0 - min_Cm0)/min_Cm0); 

% Alternative computation of Cm0 (mean value of Cm0 function) (unused)
Cm0_2 = 0;
for i = i1:1:i2-1
    Cm0_2 = Cm0_2 + (Cl_DVM(i+1)-Cl_DVM(i))*(Cm0_DVM(i+1)+Cm0_DVM(i))/2;
end
Cm0_2 = Cm0_2/(Cl_DVM(i2)-Cl_DVM(i1));






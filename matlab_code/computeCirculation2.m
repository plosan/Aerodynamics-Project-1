function [x, z, vortex, node, c, n_vec, t_vec, Gamma] = ...
    computeCirculation2(NACA, U_inf, chord, x_flap, eta, alpha, N, distribution)
% Computes geometry and circulation
%--------------------------------------------------------------------------
% Inputs:
%   - Physical data
%       - NACA          Airfoil (0000 - 9999)
%       - U_inf         Free stream velocity
%       - chord         Airfoil chord
%       - x_flap        Flap position
%       - eta           Flap deflection angle
%       - alpha         Attack angle
%   - Numerical data
%       - N             Number of nodes for panels (N-1 panels)
%       - distribution  Type of distribution (uniform or full cosine)
%--------------------------------------------------------------------------
% Outputs:
%   - x         Nodes position (x coordinate)           
%   - z         Nodes position (z coordinate)
%   - vortex    Vortex position [N-1 x 2]
%   - node      Control nodes position [N-1 x 2]
%   - c         Discretization panel chords
%   - n_vec     Normal vectors [N-1 x 2]
%   - t_vec     Tangent vectors [N-1 x 2]
%   - Gamma     Circulation [N-1 x 1]
%--------------------------------------------------------------------------

% Compute airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)
%t = mod(NACA, 100)/100;             % Thickness (percent of chord)

% Compute distribution (linear or full cosine)
if distribution == "uniform"
    x = linspace(0, 1, N);         % Linear distribution
else
    theta = linspace(0, pi, N);    
    x = (1 - cos(theta))/2;             % Cosine distribution
end

% Compute camber line coordinates (Z)
z = zeros(1, N);          % Camber line
i = 1;
while x(i) < p
    z(i) = f/p^2*(2*p*x(i) - x(i)^2);
    i = i+1;
end
z(i:end) = f/(1-p)^2*(1 - 2*p + 2*p*x(i:end) - x(i:end).^2);
x = x*chord;
z = z*chord;

% Find closest point to flap
min_dif = intmax;
ih = 1;
for i = 1:N
    dif = abs(x_flap - x(i));
    if dif < min_dif
        min_dif = dif;
        ih = i;
    end
end

xh = x(ih);
zh = z(ih);

% Rotate flap
x(ih:end) = x(ih:end) - x(ih);  % Move to origin (x)
z(ih:end) = z(ih:end) - z(ih);  % Move to origin (z)

M = [cos(eta) sin(eta); -sin(eta) cos(eta)];    % Rotation matrix
X = M*[x(ih:end); z(ih:end)];                   % Rotation
                                   
x(ih:end) = X(1,:);
z(ih:end) = X(2,:);

x(ih:end) = x(ih:end) + xh;     % Move to xh
z(ih:end) = z(ih:end) + zh;     % Move to zh


% Compute geometric properties of discretization
c = zeros(1, N-1);
t_vec = zeros(N-1, 2);
n_vec = t_vec;
vortex = t_vec;
node = t_vec;

for i = 1:N-1
    delta_x = x(i+1) - x(i);
    delta_z = z(i+1) - z(i);
    c(i) = sqrt(delta_x^2 + delta_z^2);
    t_vec(i,:) = [delta_x delta_z]/c(i);
    n_vec(i,:) = [-delta_z delta_x]/c(i);
    vortex(i,:) = [x(i) z(i)] + t_vec(i,:)*c(i)/4;
    node(i,:) = [x(i) z(i)] + 3*t_vec(i,:)*c(i)/4;
end

% Compute influence coefficients and RHS
A = zeros(N-1,N-1);
RHS = zeros(N-1,1);
for i = 1:N-1
    for j = 1:N-1
        r_sq = (node(i,1)-vortex(j,1))^2 + (node(i,2)-vortex(j,2))^2;
        u = (node(i,2)-vortex(j,2))/(2*pi*r_sq);
        w = -(node(i,1)-vortex(j,1))/(2*pi*r_sq);
        A(i,j) = n_vec(i,1)*u + n_vec(i,2)*w;
    end
    RHS(i) = -U_inf*(n_vec(i,1)*cos(alpha) + n_vec(i,2)*sin(alpha));
end

% Solve linear system
Gamma = linsolve(A, RHS);

end
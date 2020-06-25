function [Gamma] = computeCirculation(U_inf, alpha, vortex, node, n_vec, N)
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

% Compute influence coefficients and RHS
A = zeros(N, N);
RHS = zeros(N,1);
for i = 1:N
    for j = 1:N
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
function [Gamma] = computeCirculation(U_inf, alpha, vortex, node, n_vec, N)
%--------------------------------------------------------------------------
% Inputs:
% Physical data
%       - U_inf     Free stream velocity
%       - alpha     Attack angle
%       - vortex    Vortex position [N x 2]
%       - node      Control nodes position [N x 2]
%       - n_vec     Normal vectors [N x 2]
% Numerical data
%       - N        Number of panels (N+1 nodes)
%--------------------------------------------------------------------------
% Outputs:
%   - Gamma     Circulation [N x 1]
%--------------------------------------------------------------------------

% Compute influence coefficients and RHS
A = zeros(N, N);    % Matrix of influence coefficients
RHS = zeros(N,1);   % Right hand side vector for DVM
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
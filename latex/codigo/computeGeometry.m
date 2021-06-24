function [x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N, distribution)
% Compute airfoil geometry (camber line, vortex nodes, control nodes,
% panel chord, normal vectors)

% Compute distribution (uniform or full cosine)
if distribution == "uniform"    % Uniform distribution
    x = linspace(0, 1, N+1);   
else % Full cosine distribution
    theta = linspace(0, pi, N+1);    
    x = (1 - cos(theta))/2; 
end

% Compute camber line coordinates (Z)
z = zeros(1, N+1);          % Camber line
i = 1;
while x(i) < p % Prevent division by 0, in case of p = 0
    z(i) = f/p^2*(2*p*x(i) - x(i)^2);
    i = i+1;
end
z(i:end) = f/(1-p)^2*(1 - 2*p + 2*p*x(i:end) - x(i:end).^2);
x = x*chord;
z = z*chord;

% Find closest point to flap
if x_flap ~= 0 % in case there is a trailing edge flap
    min_dif = intmax;
    ih = 1;
    for i = 1:N+1
        dif = abs(x_flap - x(i));
        if dif < min_dif
            min_dif = dif;
            ih = i;
        end
    end
    % Flap initial node coordinates
    xh = x(ih);
    zh = z(ih);
    % Rotate flap
    x(ih:end) = x(ih:end) - x(ih);  % Move to origin (x)
    z(ih:end) = z(ih:end) - z(ih);  % Move to origin (z)
    M = [cos(eta) sin(eta); -sin(eta) cos(eta)];    % Rotation matrix
    X = M*[x(ih:end); z(ih:end)];                   % Rotation
    x(ih:end) = X(1,:); % New flap coordinates (x)
    z(ih:end) = X(2,:); % New flap coordinates (z)
    x(ih:end) = x(ih:end) + xh;     % Move to xh
    z(ih:end) = z(ih:end) + zh;     % Move to zh
end

% Compute geometric properties of discretization
c = zeros(1, N);        % Panel's chord
t_vec = zeros(N, 2);    % Tangent vectors (unused)
n_vec = t_vec;          % Normal vectors
vortex = t_vec;         % Vortex position
node = t_vec;           % Control nodes (solid body condition)
for i = 1:N
    delta_x = x(i+1) - x(i);
    delta_z = z(i+1) - z(i);
    c(i) = sqrt(delta_x^2 + delta_z^2);
    t_vec(i,:) = [delta_x delta_z]/c(i);
    n_vec(i,:) = [-delta_z delta_x]/c(i);
    vortex(i,:) = [x(i) z(i)] + t_vec(i,:)*c(i)/4;
    node(i,:) = [x(i) z(i)] + 3*t_vec(i,:)*c(i)/4;
end

end
function [x, z, vortex, node, c, n_vec, t_vec] = ...
    computeGeometry(f, p, chord, x_flap, eta, N, distribution)

% Compute distribution (linear or full cosine)
if distribution == "uniform"
    x = linspace(0, 1, N+1);         % Linear distribution
else
    theta = linspace(0, pi, N+1);    
    x = (1 - cos(theta))/2;             % Cosine distribution
end

% Compute camber line coordinates (Z)
z = zeros(1, N+1);          % Camber line
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
for i = 1:N+1
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
c = zeros(1, N);
t_vec = zeros(N, 2);
n_vec = t_vec;
vortex = t_vec;
node = t_vec;

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
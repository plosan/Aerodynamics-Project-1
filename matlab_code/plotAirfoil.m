function plotAirfoil(NACA, alpha, points, distribution)

% Compute airfoil geometric parameters
f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)
t = mod(NACA, 100)/100;             % Thickness (percent of chord)

% alpha to radians
alpha = alpha*pi/180;

% Compute distribution (linear or cosine)
if distribution == "linear"
    x = linspace(0, 1, points);         % Linear distribution
else
    theta = linspace(0, pi, points);    
    x = (1 - cos(theta))/2;             % Cosine distribution
end

% Compute camber line coordinates (zc) and derivative of zc respect to x
% (tangent)
zc = zeros(1, points);          % Camber line
der_zc_x = zeros(1, points);    % (d zc)/(dx)
i = 1;
while x(i) <= p
    zc(i) = f/p^2*(2*p*x(i) - x(i)^2);
    der_zc_x(i) = 2*f/p^2*(p - x(i));
    i = i+1;
end
zc(i:end) = f/(1-p)^2*(1 - 2*p + 2*p*x(i:end) - x(i:end).^2); 
der_zc_x(i:end) = 2*f/(1-p)^2*(p - x(i:end));

% Compute thickness distribution
zt = t/0.2*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4);

% Compute sine and cosine of camber line angle
cos_theta = 1./sqrt(1 + der_zc_x.^2);
sin_theta = der_zc_x./sqrt(1 + der_zc_x.^2); 

% Compute coordinates of the upper part of the airfoil
xu = x - zt.*sin_theta;
zu = zc + zt.*cos_theta;

% Compute coordinates of the lower part of the airfoil
xl = x + zt.*sin_theta;
zl = zc - zt.*cos_theta;

% Rotate 
if alpha ~= 0
    % Rotation matrix
    M = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
    % Rotate camber line
    X = M*[x; zc];
    x = X(1,:);
    zc = X(2,:);
    % Rotate upper line
    X = M*[xu; zu];
    xu = X(1,:);
    zu = X(2,:);
    % Rotate lower line
    X = M*[xl; zl];
    xl = X(1,:);
    zl = X(2,:);    
end

% Plot Airfoil
str = sprintf("\\textbf{NACA %d, $\\alpha = %.1f ^\\circ$}", NACA, alpha*180/pi);
figure;
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title(str);
plot(x, zc, 'r');
plot(xu, zu, 'b');
plot(xl, zl, 'b');
axis equal;
if alpha == 0
    xlim([0 1]);
    ylim([-0.5 0.5]);
end
xlabel("$x / c$");
ylabel("$z / c$");
grid on;
box on;
set(gcf, 'units', 'centimeters', 'position', [1,1,18,15]);
hold off;

end
clear;
close all;
clc;

NACA = 4412;

f = floor(NACA/1000)/100;           % Maximum camber (percent of chord)
p = mod(floor(NACA/100), 10)/10;    % Maximum camber position (tenths of chord)
t = mod(NACA, 100)/100;             % Thickness (percent of chord)
alpha = 4*pi/180;

[A0, A1, A2] = computeACoefficients(f, p, alpha);


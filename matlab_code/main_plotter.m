clear;
close all;
clc;

points = 1e3;

NACA = 0000;
plotAirfoil(NACA, 0, points, "cosine")

NACA = 2206;
plotAirfoil(NACA, 0, points, "cosine")





clear;
close all;
clc;

%% 1. DATA INPUT

[NACA, f, p, t, alpha, xh, eta, N] = userInput();

% Print input data
fprintf("NACA Airfoil: %d\n", NACA);
fprintf("Max camber (%% chord): %.2f\n", f);
fprintf("Max camber position (%% chord): %.2f\n", p);
fprintf("Relative thickness (%% chord): %.2f\n", f);
fprintf("Angle of attack: %.2f\n", alpha);
fprintf("Flap hinge location: %.2f\n", xh);
fprintf("Flap deflection angle: %.2f\n", eta);
fprintf("Number of stages: %d\n", N);

fprintf("\n\n");
fprintf("%35s: %d\n", "NACA Airfoil", NACA);
fprintf("%35s: %.2f\n", "Max camber (%% chord)", f);
fprintf("%35s: %.2f\n", "Max camber position (1/10 chord)", p);
fprintf("%35s: %.2f\n", "Relative thickness (% chord)", f);
fprintf("%35s: %.2f\n", "Angle of attack (degrees)", alpha);
fprintf("%35s: %.2f\n", "Flap hinge location (% chord)", xh);
fprintf("%35s: %.2f\n", "Flap deflection angle (degrees)", eta);
fprintf("%35s: %d\n", "Number of stages", N);




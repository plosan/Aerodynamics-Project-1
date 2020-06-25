function plotAirfoilDVM(NACA, x, z, vortex, node, n_vec)

% Plot Airfoil
figure();
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title(sprintf("\\textbf{NACA %d}", NACA));
plot(x, z, 'r', 'LineWidth', 1);
scatter(x, z, 20, 'r', 'filled');
scatter(vortex(:,1), vortex(:,2), 20, 'g', 'filled');
scatter(node(:,1), node(:,2), 20, 'b', 'filled');
quiver(node(:,1), node(:,2), n_vec(:,1), n_vec(:,2), 'LineWidth', 1);
axis equal;
xlim([-0.1 1.1]);
ylim([-0.1 0.5]);
xlabel("$x / c$");
ylabel("$z / c$");
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.1f'));
grid on;
box on;
set(gcf, 'units', 'centimeters', 'position', [18,1,18,10]);
legend("L\'inea media", "Puntos de discretizaci\'on", "Nodos de v\'ortice", ...
    "Nodos de control", "Vectores normales");
hold off;

end
function plotDVMCoefficients(alpha, Cl, Cm0, i1, i2)

alpha = alpha*180/pi;

figure();
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Coeficiente de sustentaci\'on $\left( C_l \right)$}");
plot(alpha, Cl, 'b');
scatter(alpha(i1), Cl(i1), 20, 'r', 'filled');
scatter(alpha(i2), Cl(i2), 20, 'r', 'filled');
% xlabel("\'Angulo de ataque $\left( ^\circ \right)$");
% ylabel("Coeficiente de sustentaci\'on");
xlabel("$\alpha \ \left( \mathrm{deg} \right)$");
ylabel("$C_l$");
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.1f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.1f'));
set(gca, 'XAxisLocation', 'origin');
set(gca, 'YAxisLocation', 'origin');
set(gcf, 'units', 'centimeters', 'position', [0,1,18,15]);
legend("Coeficiente de sustentaci\'on $\left( C_l \right)$", ...
    "L\'imites r\'egimen lineal");
grid on;
grid minor;
box on;
hold off;

figure();
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Coeficiente de momento libre $\left( C_{m0} \right)$}");
plot(Cl, Cm0, 'b');
scatter(Cl(i1), Cm0(i1), 20, 'r', 'filled');
scatter(Cl(i2), Cm0(i2), 20, 'r', 'filled');
xlabel("Coeficiente de sustentaci\'on");
ylabel("Coeficiente de momento libre");
ylim([-0.1 0]);
set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.2f'));
set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
set(gcf, 'units', 'centimeters', 'position', [18,1,18,15]);
legend("Coeficiente de momento libre $\left( C_{m0} \right)$", ...
    "L\'imites r\'egimen lineal");
grid on;
grid minor;
box on;
hold off;



end
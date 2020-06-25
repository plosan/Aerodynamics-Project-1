function plotError(N, Cl_TAT, Cl_DVM, Cm_TAT, Cm_DVM)

%% LIFT COEFFICIENT PLOT
figure(1);
hold on;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
title("\textbf{Coeficiente de sustetaci\'on $\left( C_l \right)$ seg\'un TAT, DVM y Error relativo}");
% Left axis
    yyaxis left;
    plot([0 max(N)], [Cl_TAT, Cl_TAT], '-k');
    plot(N, Cl_DVM, '-b');
    text(20, Cl_TAT+0.01, sprintf("\\textbf{$C_{l,\\mathrm{TAT}} = %.4f$}", Cl_TAT));
    xlabel("N\'umero de paneles");
    ylabel("Coeficiente de sustentaci\'on");
    set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.0f'));
    set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
    set(gca,'ycolor','k'); 
% Right axis
    yyaxis right;
    err_Cl = 100*abs((Cl_TAT-Cl_DVM)/Cl_TAT);
    plot(N, err_Cl, '-r');
    ylabel("Error relativo $\left( \% \right)$");
    set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.0f'));
    set(gca,'ycolor','k');
set(gcf, 'units', 'centimeters', 'position', [0,1,18,11]);
legend("$C_l$ TAT", "$C_l$ DVM", "Error relativo $C_l$", "Location", "Southeast");
grid on;
grid minor;
box on;
hold off;

%% MOMENT COEFFICIENT PLOT
figure(2);
hold on;
title("\textbf{Coeficiente de momentos LE $\left( C_{m,\mathrm{LE}} \right)$ seg\'un TAT, DVM y Error relativo}");
% Left axis
    yyaxis left;
    plot([0 max(N)], [Cm_TAT Cm_TAT], '-k');
    plot(N, Cm_DVM, '-b');
    text(20, Cm_TAT-0.005, sprintf("\\textbf{$C_{m,\\mathrm{TAT}} = %.4f$}", Cm_TAT));
    xlabel("N\'umero de paneles");
    ylabel("Coeficiente de momento LE");
    ylim([-0.24 -0.08]);
    set(gca, 'xticklabel', num2str(get(gca,'xtick')', '%.0f'));
    set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.2f'));
    set(gca,'ycolor', 'k');
% Right axis
    yyaxis right;
    err_Cm = 100*(abs(Cm_DVM-Cm_TAT)/Cm_TAT);
    plot(N, err_Cm, 'r');
    ylabel("Error relativo $\left( \% \right)$");
    set(gca, 'yticklabel', num2str(get(gca,'ytick')', '%.0f'));
    set(gca,'ycolor','k');
set(gcf, 'units', 'centimeters', 'position', [18,1,18,11]);
legend("$C_{m, \mathrm{LE}}$ TAT", "$C_{m, \mathrm{LE}}$ DVM", "Error relativo $C_{m, \mathrm{LE}}$", "Location", "Northeast");
grid on;
grid minor;
box on;
hold off;

end
% plotting relative permeability and capillary pressure curves
figure(1)
kr = load('Materials/krCurveSand_200.dat');
plot(kr(1:89,1),kr(1:89,2), 'LineWidth',1.5,'Color','b')
xlabel('Capillary Pressure p (kPa)')
ylabel('Relative permeability k_r')
xlim([0 10])

set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 16);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 10)
stmp = strcat('Images\', 'Rel_perm', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

figure(2)
pc = load('Materials/pcCurveSand_200.dat');
plot(pc(:,1),pc(:,2), 'LineWidth',1.5,'Color','r')
xlabel('Capillary pressure p (kPa)')
xlim([0 40])
ylabel('Effective Saturation S^*_w')

set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 16);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 10)
stmp = strcat('Images\', 'Effective_sat', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

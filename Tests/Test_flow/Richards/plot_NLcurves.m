% plotting relative permeability and capillary pressure curves
figure(1)
kr = load('Materials/krCurveSand_200.dat');
plot(kr(1:89,1),kr(1:89,2), 'LineWidth',1,'Color','b')
xlabel('Pressure p (kPa)')
ylabel('Relative permeability k_r')
xlim([0 10])

set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif')


figure(2)
pc = load('Materials/pcCurveSand_200.dat');
plot(pc(:,1),pc(:,2), 'LineWidth',1,'Color','r')
xlabel('Capillary pressure p (kPa)')
xlim([0 40])
ylabel('Effective Saturation S^*_w')

set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif')
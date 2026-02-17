function plotPQ(p, q, A, B, druckerPrager)

	% Given material parameters
	cohes_DZ1 = druckerPrager.DamageZone1.cohesion_mean;
	phi_DZ1 = druckerPrager.DamageZone1.friction_angle_mean;

	cohes_core = druckerPrager.Core.cohesion_mean;
	phi_core = druckerPrager.Core.friction_angle_mean;

	cohes_DZ2 = druckerPrager.DamageZone2.cohesion_mean;
	phi_DZ2 = druckerPrager.DamageZone2.friction_angle_mean;

	% --- Damage Zone 1 ---
	phi_DZ1 = deg2rad(phi_DZ1);
	A_DZ1 = -3.0*tan(phi_DZ1) / sqrt(9.0 + 12.0*tan(phi_DZ1)^2);
	B_DZ1 = 3.0*cohes_DZ1 / sqrt(9.0 + 12.0*tan(phi_DZ1)^2);

	% --- Core ---
	phi_core = deg2rad(phi_core);
	A_core = -3.0*tan(phi_core) / sqrt(9.0 + 12.0*tan(phi_core)^2);
	B_core = 3.0*cohes_core / sqrt(9.0 + 12.0*tan(phi_core)^2);

	% --- Damage Zone 2 ---
	phi_DZ2 = deg2rad(phi_DZ2);
	A_DZ2 = -3.0*tan(phi_DZ2) / sqrt(9.0 + 12.0*tan(phi_DZ2)^2);
	B_DZ2 = 3.0*cohes_DZ2 / sqrt(9.0 + 12.0*tan(phi_DZ2)^2);

	% Plot
	figure;
	hold on;

	plot(p, A_DZ1*p + B_DZ1, 'LineWidth', 2, 'DisplayName', ...
		sprintf('damage zone 1: q = %.3f p + %.3f', A_DZ1, B_DZ1));

	plot(p, A_core*p + B_core, 'LineWidth', 2, 'DisplayName', ...
		sprintf('core: q = %.3f p + %.3f', A_core, B_core));

	plot(p, A_DZ2*p + B_DZ2, 'LineWidth', 2, 'DisplayName', ...
		sprintf('damage zone 2: q = %.3f p + %.3f', A_DZ2, B_DZ2));

	plot(p, q, 'o--', 'LineWidth', 2, 'DisplayName', ...
		sprintf('average: q = %.3f p + %.3f', A, B));

	xlabel('p [MPa]');
	ylabel('q [MPa]');
	grid on;
	legend('Location','northeast');

	hold off;

	% Save figure
	set(gcf,'PaperPositionMode','auto');
	print(gcf, 'pq_plot.pdf', '-dpdf', '-bestfit');

end
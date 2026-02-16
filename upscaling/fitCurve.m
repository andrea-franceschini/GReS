function fitCurve(avgPs, avgQs)

    % Linear fit q = A*p + B
    p = polyfit(avgPs, avgQs, 1);
    A = p(1);
    B = p(2);

    % Residual
    residual = norm(polyval(p,avgPs) - avgQs)^2;

    tan_phi = -3.0*A / sqrt(9.0 - 12.0*A^2);
    phi = rad2deg(atan(tan_phi));
    cohes = B * sqrt(9.0 + 12.0*tan_phi^2) / 3.0;

    if residual < 1e-2
        fprintf('Linear model\n');
        fprintf('q = A*p + B with A %15.6e and B %15.6e\n', A, B);
    else
        fprintf('Non linear model. Residual %15.6e\n', residual);
        fprintf('Best linear fit: q = A*p + B with A %15.6e and B %15.6e\n', A, B);
    end

    fprintf('Equivalent cohesion and friction angle are %15.6e and %15.6e\n', ...
            cohes, phi);

end

clear
close all
clc

scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);

fname = fullfile('Input','pressurizedCrack.xml');

%% ========================================================================
%  SNEDDON TEST PARAMETERS
%  ========================================================================

% Domain dimensions
Lx = 100.0;
Ly = 1.0;
Lz = 100.0;

% Base grid used by Cusini: 15 x 15.
% Uniform 3x3 refinements:
%
%   15 -> 45 -> 135 -> 405
%
% Add 405 if desired.
nList = [15, 45, 135, 405];

% The fracture spans exactly 7 cells of the 15x15 base grid.
% DO NOT use 46.67: the exact value keeps the fracture tips on grid
% boundaries at every refinement level.
fractureLength = 7*(Lx/15);
a = fractureLength/2;

% Mechanical parameters
%
% IMPORTANT: these values MUST be the same as those used in the XML
% material definition.
nu = 0.25;
E  = 25000;

% Prescribed fracture traction.
% Negative = opening pressure with the current sign convention.
p = -2.0;

% Store convergence results
nRuns = numel(nList);

hVec      = zeros(nRuns,1);
errAbsVec = zeros(nRuns,1);
errRelVec = zeros(nRuns,1);

%% ========================================================================
%  REFINEMENT LOOP
%  ========================================================================

for iRun = 1:nRuns

    nCell = nList(iRun);

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf(' SNEDDON TEST: %d x %d grid\n',nCell,nCell);
    fprintf('============================================================\n');

    % Re-read input for a clean problem at every refinement
    params = readInput(fname);

    simparams = SimulationParameters(params.SimulationParameters);

    %% --------------------------------------------------------------------
    %  Grid
    %  --------------------------------------------------------------------

    grid = structuredMesh( ...
        nCell,1,nCell, ...
        [-0.5*Lx,0.5*Lx], ...
        [-0.5*Ly,0.5*Ly], ...
        [-0.5*Lz,0.5*Lz]);

    h = Lx/nCell;
    hVec(iRun) = h;

    %% --------------------------------------------------------------------
    %  Analytical Sneddon displacement on external boundary
    %
    %  Horizontal fracture:
    %
    %             z
    %             ^
    %             |
    %      --------+---------> x
    %            -a   a
    %
    %  fracture: z = 0,  -a <= x <= a
    %
    %  --------------------------------------------------------------------

    c = grid.coordinates;

    Xall = c(:,1);
    Zall = c(:,3);

    Xhalf = Lx/2;
    Zhalf = Lz/2;

    tol = 1e-10*max(Lx,Lz);

    boundNodes = find( ...
          abs(Xall - Xhalf) < tol ...
        | abs(Xall + Xhalf) < tol ...
        | abs(Zall - Zhalf) < tol ...
        | abs(Zall + Zhalf) < tol);

    X = Xall(boundNodes);
    Z = Zall(boundNodes);

    % Distances from the two fracture tips
    r1 = sqrt((X + a).^2 + Z.^2);
    r2 = sqrt((X - a).^2 + Z.^2);

    % Polar angles with respect to the tips
    theta1 = atan2(Z,X + a);
    theta2 = atan2(Z,X - a);

    % Polar coordinates with respect to the fracture centre
    r     = sqrt(X.^2 + Z.^2);
    theta = atan2(Z,X);

    sqrtR1R2 = sqrt(r1.*r2);
    theta12  = 0.5*(theta1 + theta2);

    % Exact infinite-medium displacement
    fac = p*(1 + nu)/E;

    ux = -fac .* ( ...
        (1 - 2*nu) .* ...
        (sqrtR1R2.*cos(theta12) - r.*cos(theta)) ...
        - ...
        (r.^2./sqrtR1R2).*sin(theta).* ...
        sin(theta - theta12) );

    uz = -fac .* ( ...
        2*(1 - nu) .* ...
        (sqrtR1R2.*sin(theta12) - r.*sin(theta)) ...
        + ...
        r.*sin(theta).* ...
        (1 - (r./sqrtR1R2).* ...
        cos(theta - theta12)) );

    %% --------------------------------------------------------------------
    %  Material and output
    %  --------------------------------------------------------------------

    mat = Materials(params.Materials);

    outputName = sprintf('Output/Sneddon_%d',nCell);

    printUtils = OutState( ...
        "outputFile",outputName, ...
        "printTimes",1, ...
        "matFileName",outputName);

    %% --------------------------------------------------------------------
    %  Boundary conditions
    %  --------------------------------------------------------------------

    bc = Boundaries(grid);

    % Exact x displacement on the complete x-z external boundary
    bc.addBC( ...
        'name',"x_fix", ...
        'type',"dirichlet", ...
        'field',"node", ...
        'variable',"displacements", ...
        'entityListType',"bcList", ...
        'entityList',boundNodes, ...
        'components',"x");

    bc.addBCEvent( ...
        "x_fix", ...
        'time',0.0, ...
        'value',ux);

    % Exact z displacement on the complete x-z external boundary
    bc.addBC( ...
        'name',"z_fix", ...
        'type',"dirichlet", ...
        'field',"node", ...
        'variable',"displacements", ...
        'entityListType',"bcList", ...
        'entityList',boundNodes, ...
        'components',"z");

    bc.addBCEvent( ...
        "z_fix", ...
        'time',0.0, ...
        'value',uz);

    % Plane-strain constraint in the thin y direction
    bc.addBC( ...
        'name',"y_fix", ...
        'type',"dirichlet", ...
        'field',"surface", ...
        'variable',"displacements", ...
        'entityListType',"tag", ...
        'entityList',[3,4], ...
        'components',"y");

    bc.addBCEvent( ...
        "y_fix", ...
        'time',0.0, ...
        'value',0.0);

    %% --------------------------------------------------------------------
    %  Domain
    %  --------------------------------------------------------------------

    domain = Discretizer( ...
        'Boundaries',bc, ...
        'Materials',mat, ...
        'Grid',grid);

    domain.addPhysicsSolvers(params.Solver);

    %% --------------------------------------------------------------------
    %  Prescribed fracture pressure
    %  --------------------------------------------------------------------

    efem = getPhysicsSolver(domain,"EFEMaugmented");

    % Local fracture component 1 = normal traction
    efem.bcTraction(1:3:end) = p;

    %% --------------------------------------------------------------------
    %  Solve
    %  --------------------------------------------------------------------

    solver = NonLinearImplicit( ...
        'simulationparameters',simparams, ...
        'domains',domain, ...
        'output',printUtils);

    solver.simulationLoop();

    %% ====================================================================
    %  POSTPROCESSING
    %  ====================================================================

    efem = getPhysicsSolver(domain,"EFEMaugmented");

    f = efem.fractureMesh.surfaces;

    % ---------------------------------------------------------------------
    % Fracture geometry checks
    % ---------------------------------------------------------------------

    % This convergence test is specifically for a horizontal fracture.
    %
    % Its centres should therefore lie at z = 0.
    zFrac = f.center(:,3);

    assert(max(abs(zFrac)) < 1e-8, ...
        ['This script assumes the horizontal Sneddon problem. ', ...
         'The fracture in pressurizedCrack.xml is not horizontal/centred ', ...
         'at z = 0.']);

    % Local coordinate along the horizontal fracture
    xi = f.center(:,1);

    % Numerical normal opening
    g = getState(efem,"fractureJump");
    gn = g(1:3:end);

    assert(numel(xi) == numel(gn), ...
        'Mismatch between fracture elements and fracture-jump DOFs.');

    % Sort fracture elements from left to right
    [xi,id] = sort(xi);
    gn = gn(id);

    % ---------------------------------------------------------------------
    % Verify that the fracture discretization is the expected one
    % ---------------------------------------------------------------------

    % For this benchmark:
    %
    %   base grid: 15 cells
    %   fracture:  7 base cells
    %
    % Therefore the expected number of fracture elements is
    %
    %   7 * (nCell/15)
    %
    expectedNF = round(7*nCell/15);

    if numel(xi) ~= expectedNF
        warning( ...
            ['Expected %d fracture elements but found %d.\n', ...
             'Check the fracture length in pressurizedCrack.xml.\n', ...
             'For this benchmark it should be %.15g m.'], ...
            expectedNF,numel(xi),fractureLength);
    end

    % ---------------------------------------------------------------------
    % Exact Sneddon aperture
    % ---------------------------------------------------------------------

    % Plane-strain aperture:
    %
    %   gn(x) = C*sqrt(a^2-x^2)
    %
    % with the current negative-pressure convention.
    C = -4*(1 - nu^2)*p/E;

    gnAnal = C*sqrt(max(0,a^2 - xi.^2));

    %% --------------------------------------------------------------------
    %  TRUE L2 ERROR ALONG THE FRACTURE
    %
    %  Do not simply use
    %
    %     area.*(gn - gnAnal(center)).^2
    %
    %  because the exact solution has a square-root tip singularity.
    %  Midpoint quadrature close to the tips can contaminate the measured
    %  convergence rate.
    %
    %  Since gn is constant on each EFEM(0) fracture element, integrate
    %  the error analytically over each segment.
    %  --------------------------------------------------------------------

    % Segment boundaries.
    %
    % For a horizontal fracture terminating exactly at grid boundaries,
    % the midpoint between two consecutive fracture centres is the common
    % segment boundary.
    xEdge = [-a; 0.5*(xi(1:end-1) + xi(2:end)); a];

    % Make sure all reconstructed intervals are valid
    assert(all(diff(xEdge) > 0), ...
        'Invalid ordering of fracture element centres.');

    % Primitive of sqrt(a^2-x^2)
    %
    % Integral sqrt(a^2-x^2) dx
    %
    % = 1/2 [x sqrt(a^2-x^2) + a^2 asin(x/a)]
    Fsqrt = @(x) 0.5 .* ( ...
        x .* sqrt(max(0,a^2 - x.^2)) ...
        + a^2 .* asin(max(-1,min(1,x./a))) );

    % Primitive of (a^2-x^2)
    Fsquare = @(x) a^2.*x - x.^3./3;

    err2 = 0.0;

    for k = 1:numel(gn)

        xL = xEdge(k);
        xR = xEdge(k+1);

        I0 = xR - xL;

        % Integral sqrt(a^2-x^2)
        I1 = Fsqrt(xR) - Fsqrt(xL);

        % Integral (a^2-x^2)
        I2 = Fsquare(xR) - Fsquare(xL);

        % Integral:
        %
        %   [gn_h - C sqrt(a^2-x^2)]^2 dx
        %
        err2 = err2 ...
            + gn(k)^2*I0 ...
            - 2*gn(k)*C*I1 ...
            + C^2*I2;
    end

    % Eliminate possible tiny negative roundoff
    err2 = max(err2,0);

    errAbs = sqrt(err2);

    % Exact analytical L2 norm:
    %
    % Integral_{-a}^a C^2(a^2-x^2) dx
    % = C^2 * 4*a^3/3
    exactNorm2 = C^2*(4*a^3/3);

    errRel = errAbs/sqrt(exactNorm2);

    errAbsVec(iRun) = errAbs;
    errRelVec(iRun) = errRel;

    fprintf('h                  = %1.8e\n',h);
    fprintf('fracture elements  = %d\n',numel(gn));
    fprintf('absolute L2 error  = %1.8e\n',errAbs);
    fprintf('relative L2 error  = %1.8e\n',errRel);

    %% --------------------------------------------------------------------
    %  Aperture plot
    %  --------------------------------------------------------------------

    xiAnal = linspace(-a,a,1000)';
    gnAnalPlot = C*sqrt(max(0,a^2 - xiAnal.^2));

    fig = figure(iRun);

    plot( ...
        xiAnal,gnAnalPlot, ...
        'r-', ...
        'LineWidth',1.5);

    hold on

    plot( ...
        xi,gn, ...
        'ko', ...
        'MarkerSize',3.5, ...
        'MarkerFaceColor','k');

    xlim([-1.1*a,1.1*a]);

    xlabel( ...
        '$\xi$', ...
        'Interpreter','latex', ...
        'FontSize',14);

    ylabel( ...
        '$g_N$', ...
        'Interpreter','latex', ...
        'FontSize',14);

    set( ...
       (gca), ...
        'TickLabelInterpreter','latex', ...
        'FontSize',14);

    grid on
    box on

    legend( ...
        'Analytical', ...
        sprintf('EFEM(0), $h = %.3g$',h), ...
        'Interpreter','latex', ...
        'Location','best');

    exportgraphics( ...
        fig, ...
        fullfile('Output',sprintf('gn_plot_%d.pdf',nCell)));
end

%% ========================================================================
%  OBSERVED CONVERGENCE ORDER
%  ========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf(' CONVERGENCE RESULTS\n');
fprintf('============================================================\n');

fprintf('\n');
fprintf('%10s %15s %15s %12s\n', ...
    'h','abs L2 error','rel L2 error','order');

fprintf('%10.4e %15.6e %15.6e %12s\n', ...
    hVec(1),errAbsVec(1),errRelVec(1),'-');

pObs = nan(nRuns-1,1);

for k = 2:nRuns

    pObs(k-1) = ...
        log(errRelVec(k-1)/errRelVec(k)) / ...
        log(hVec(k-1)/hVec(k));

    fprintf('%10.4e %15.6e %15.6e %12.6f\n', ...
        hVec(k), ...
        errAbsVec(k), ...
        errRelVec(k), ...
        pObs(k-1));
end

%% ========================================================================
%  CONVERGENCE PLOT
%  ========================================================================

figure

loglog( ...
    hVec,errRelVec, ...
    'ko-', ...
    'MarkerFaceColor','k', ...
    'LineWidth',1.5);

hold on

% First-order reference line through the finest-grid result
ref = errRelVec(end)*(hVec/hVec(end));

loglog( ...
    hVec,ref, ...
    'k--', ...
    'LineWidth',1.0);

xlabel( ...
    '$h$', ...
    'Interpreter','latex', ...
    'FontSize',14);

% ylabel( ...
%     '$\|g_{N,h}-g_N^{\mathrm{an}}\|_{L^2}/' + ...
%     '\|g_N^{\mathrm{an}}\|_{L^2}$', ...
%     'Interpreter','latex', ...
%     'FontSize',14);

legend( ...
    'EFEM(0)', ...
    '$\mathcal{O}(h)$', ...
    'Interpreter','latex', ...
    'Location','best');

set( ...
   (gca), ...
    'TickLabelInterpreter','latex', ...
    'FontSize',14);

grid on
box on

exportgraphics( ...
    gcf, ...
    fullfile('Output','Sneddon_convergence.pdf'));
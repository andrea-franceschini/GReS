function anal = validateDruckerPragerAnalytical(dpLaw, driverOut, varargin)
%VALIDATEDRUCKERPRAGERANALYTICAL Semi-analytical triaxial validation for Drucker-Prager.
%
%   anal = validateDruckerPragerAnalytical(dpLaw, driverOut)
%   anal = validateDruckerPragerAnalytical(dpLaw, driverOut, 'Name', value, ...)
%
%   Inputs
%   ------
%   dpLaw     : DruckerPrager material object.
%   driverOut : matrix returned by ConstitutiveDriver.
%               Expected default columns:
%                 1 TIME
%                 2 SIG0  axial stress
%                 3 SIG1  radial stress
%                 4 SIG2  radial stress copy
%                 5 EPS0  axial strain
%                 6 EPS1  radial strain
%                 7 EPS2  radial strain copy
%                 8 ITER
%                 9 NORM
%                10 EXITFLAG
%
%   Name-value options
%   ------------------
%   'Columns'       : struct with fields time,sig0,sig1,sig2,eps0,eps1,eps2.
%   'InitialStress' : initial isotropic stress. Default: driverOut(1,SIG0).
%   'NumSubsteps'   : analytical substeps per driver increment. Default: 20.
%   'MakePlots'     : true/false. Default: true.
%   'StressScale'   : plotting scale for stress. Default: 1e-6.
%   'StressUnit'    : plotting unit label. Default: 'MPa'.
%
%   The semi-analytical solution follows the GEOS Drucker-Prager triaxial
%   driver validation formula for mixed control: prescribed axial strain and
%   prescribed radial confining stress.

  opt.Columns = struct('time',1,'sig0',2,'sig1',3,'sig2',4, ...
                       'eps0',5,'eps1',6,'eps2',7);
  opt.InitialStress = driverOut(1,opt.Columns.sig0);
  opt.NumSubsteps = 20;
  opt.MakePlots = true;
  opt.StressScale = 1e-6;
  opt.StressUnit = 'MPa';

  opt = readInput(opt,varargin{:});

  col = opt.Columns;

  timeNum      = driverOut(:,col.time);
  axStrainNum  = driverOut(:,col.eps0);
  raStrainNum  = driverOut(:,col.eps1);
  axStressNum  = driverOut(:,col.sig0);
  raStressNum  = driverOut(:,col.sig1);

  E  = dpLaw.E;
  nu = dpLaw.nu;
  if numel(E) > 1,  E  = E(1);  end
  if numel(nu) > 1, nu = nu(1); end

  cohesion      = dpLaw.co;
  frictionAngle = dpLaw.phi;
  dilationAngle = dpLaw.psi;
  hardeningRate = dpLaw.h;

  if numel(cohesion) > 1,      cohesion = cohesion(1);      end
  if numel(frictionAngle) > 1, frictionAngle = frictionAngle(1); end
  if numel(dilationAngle) > 1, dilationAngle = dilationAngle(1); end
  if numel(hardeningRate) > 1, hardeningRate = hardeningRate(1); end

  G = E/(2*(1 + nu));
  K = E/(3*(1 - 2*nu));
  lambda = K - 2*G/3;

  phi = deg2rad(frictionAngle);
  psi = deg2rad(dilationAngle);

  a = 6*cohesion*cos(phi)/(3 - sin(phi));
  b = 6*sin(phi)/(3 - sin(phi));
  bDilation = 6*sin(psi)/(3 - sin(psi));

  plasticYoungModulus = 1/(1/E + (bDilation - 3)*(b - 3)/(9*hardeningRate));
  plasticModulusForRaStrain = 1/(1/(2*G) - (b - 3)/(2*hardeningRate));

  plasticYoungModulusUnload = 1/(1/E + (bDilation + 3)*(b + 3)/(9*hardeningRate));
  plasticModulusForRaStrainUnload = 1/(1/(2*G) + (b + 3)/(2*hardeningRate));

  tAnal = [];
  axStrainAnal = [];
  raStressAnal = [];

  for n = 1:numel(timeNum)-1
    localTime = linspace(timeNum(n), timeNum(n+1), opt.NumSubsteps + 1).';
    if n > 1
      localTime = localTime(2:end);
    end

    tAnal = [tAnal; localTime]; %#ok<AGROW>
    axStrainAnal = [axStrainAnal; interp1(timeNum, axStrainNum, localTime, 'linear')]; %#ok<AGROW>
    raStressAnal = [raStressAnal; interp1(timeNum, raStressNum, localTime, 'linear')]; %#ok<AGROW>
  end

  nAnal = numel(tAnal);
  axStressAnal = zeros(nAnal,1);
  raStrainAnal = zeros(nAnal,1);

  axStressAnal(1) = opt.InitialStress;
  raStrainAnal(1) = raStrainNum(1);

  for i = 2:nAnal
    dAxStrain = axStrainAnal(i) - axStrainAnal(i-1);
    dRaStress = raStressAnal(i) - raStressAnal(i-1);

    dRaStrain = (dRaStress - lambda*dAxStrain)/(2*lambda + 2*G);
    dAxStress = (lambda + 2*G)*dAxStrain + ...
                lambda/(lambda + G)*(dRaStress - lambda*dAxStrain);

    axStressTrial = axStressAnal(i-1) + dAxStress;
    raStrainTrial = raStrainAnal(i-1) + dRaStrain;
    raStress = raStressAnal(i);

    pTrial = (axStressTrial + 2*raStress)/3;
    qTrial = -(axStressTrial - raStress);

    if qTrial >= 0
      yieldFunction = qTrial + b*pTrial - a;

      if yieldFunction >= 0
        dAxStress = dAxStrain*plasticYoungModulus;
        dRaStrain = dAxStrain - dAxStress/plasticModulusForRaStrain;

        axStressTrial = axStressAnal(i-1) + dAxStress;
        raStrainTrial = raStrainAnal(i-1) + dRaStrain;

        da = (b - 3)/3*dAxStress;
        a = a + da;
      end
    else
      yieldFunction = -qTrial + b*pTrial - a;

      if yieldFunction >= 0
        dAxStress = dAxStrain*plasticYoungModulusUnload;
        dRaStrain = dAxStrain - dAxStress/plasticModulusForRaStrainUnload;

        axStressTrial = axStressAnal(i-1) + dAxStress;
        raStrainTrial = raStrainAnal(i-1) + dRaStrain;

        da = (b + 3)/3*dAxStress;
        a = a + da;
      end
    end

    axStressAnal(i) = axStressTrial;
    raStrainAnal(i) = raStrainTrial;
  end

  pAnal = (axStressAnal + 2*raStressAnal)/3;
  qAnal = -(axStressAnal - raStressAnal);
  volStrainAnal = axStrainAnal + 2*raStrainAnal;

  pNum = (axStressNum + 2*raStressNum)/3;
  qNum = -(axStressNum - raStressNum);
  volStrainNum = axStrainNum + 2*raStrainNum;

  axStressAtNum = interp1(tAnal, axStressAnal, timeNum, 'linear');
  raStrainAtNum = interp1(tAnal, raStrainAnal, timeNum, 'linear');
  qAtNum = interp1(tAnal, qAnal, timeNum, 'linear');
  volStrainAtNum = interp1(tAnal, volStrainAnal, timeNum, 'linear');

  anal = struct();
  anal.time = tAnal;
  anal.axialStrain = axStrainAnal;
  anal.radialStrain = raStrainAnal;
  anal.axialStress = axStressAnal;
  anal.radialStress = raStressAnal;
  anal.p = pAnal;
  anal.q = qAnal;
  anal.volumetricStrain = volStrainAnal;
  anal.error.axialStress = axStressNum - axStressAtNum;
  anal.error.radialStrain = raStrainNum - raStrainAtNum;
  anal.error.q = qNum - qAtNum;
  anal.error.volumetricStrain = volStrainNum - volStrainAtNum;
  anal.error.relativeAxialStressL2 = norm(anal.error.axialStress)/max(norm(axStressAtNum),eps);
  anal.error.relativeRadialStrainL2 = norm(anal.error.radialStrain)/max(norm(raStrainAtNum),eps);
  anal.error.relativeQL2 = norm(anal.error.q)/max(norm(qAtNum),eps);
  anal.error.relativeVolumetricStrainL2 = norm(anal.error.volumetricStrain)/max(norm(volStrainAtNum),eps);

  fprintf('Drucker-Prager semi-analytical validation\n');
  fprintf('  rel. L2 axial stress      : %.4e\n', anal.error.relativeAxialStressL2);
  fprintf('  rel. L2 radial strain     : %.4e\n', anal.error.relativeRadialStrainL2);
  fprintf('  rel. L2 q                 : %.4e\n', anal.error.relativeQL2);
  fprintf('  rel. L2 volumetric strain : %.4e\n', anal.error.relativeVolumetricStrainL2);

  if opt.MakePlots
    figure('Name','Drucker-Prager triaxial validation');
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(-axStrainNum*100, qNum*opt.StressScale, 'o');
    hold on;
    plot(-raStrainNum*100, qNum*opt.StressScale, 'o');
    plot(-axStrainAnal*100, qAnal*opt.StressScale, '-','LineWidth',1.8);
    plot(-raStrainAnal*100, qAnal*opt.StressScale, '-','LineWidth',1.8);
    grid on;
    xlabel('Strain (%)');
    ylabel(['Deviatoric stress (' opt.StressUnit ')']);
    legend('num axial','num radial','anal axial','anal radial','Location','best');

    nexttile;
    plot(-axStrainNum*100, -volStrainNum*100, 'o');
    hold on;
    plot(-axStrainAnal*100, -volStrainAnal*100, '-','LineWidth',1.8);
    grid on;
    xlabel('Axial strain (%)');
    ylabel('Volumetric strain (%)');
    legend('driver','semi-analytical','Location','best');

    nexttile;
    plot(-pNum*opt.StressScale, qNum*opt.StressScale, 'o');
    hold on;
    plot(-pAnal*opt.StressScale, qAnal*opt.StressScale, '-','LineWidth',1.8);
    grid on;
    xlabel(['Mean stress (' opt.StressUnit ')']);
    ylabel(['Deviatoric stress (' opt.StressUnit ')']);
    legend('driver','semi-analytical','Location','best');
  end
end
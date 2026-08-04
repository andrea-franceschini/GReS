function testSoftSoilCreep()
% TESTSOFTSOILCREEP Unit test and sanity check for SoftSoilCreep class.
%
% Verifies initialization, stress integration, preconsolidation update,
% and tangent operator symmetry and positive-definiteness under GReS
% tension-positive / compression-negative stress conventions.

fprintf('Running SoftSoilCreep unit test...\n');

% 1. Define valid material parameter set (Isotton et al., 2019)
params = struct( ...
  'poissonRatio', 0.30, ...
  'modifiedCompressionIndex', 0.0550, ...
  'modifiedSwellingIndex', 0.0071, ...
  'creepIndex', 0.0015, ...
  'referenceTime', 86400.0, ...       % 1 day [s]
  'cohesion', 0.0, ...
  'frictionAngle', 30.0, ...          % [degrees]
  'K0NC', 0.43, ...
  'initialOCR', 1.5, ...
  'maxIterations', 50, ...
  'absoluteTolerance', 1e-12, ...
  'relativeTolerance', 1e-8);

% 2. Instantiate SoftSoilCreep object
mat = SoftSoilCreep(params);

% 3. Set up initial stress state (GReS convention: negative compression)
%    sigma = [sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_zx]
sigma0 = -[0.43, 0.43, 1.00, 0.0, 0.0, 0.0] * 1e6; % [Pa]

% Compressive strain increment in z-direction
deps = -[0.0, 0.0, 1e-5, 0.0, 0.0, 0.0];
dt = 3600.0; % 1 hour time step [s]

% 4. Test Status Initialization
status0 = mat.initializeStatus(sigma0);

assert(size(status0, 1) == 1 && size(status0, 2) == 2, ...
  'SoftSoilCreepTest:InvalidStatusShape', 'Initial status shape must be [nGauss x 2].');
assert(status0(1,1) > 0, ...
  'SoftSoilCreepTest:InvalidPreconsolidation', 'Initial preconsolidation pressure p_c,r must be positive.');

% 5. Test Local Stress Integration and Tangent Matrix Generation
cellID = 1;
[D, sigma1, status1] = mat.getStiffnessMatrix(sigma0, deps, dt, status0, cellID);

% 6. Numerical Stability and Finite Output Checks
assert(all(isfinite(sigma1(:))), ...
  'SoftSoilCreepTest:NonFiniteStress', 'Updated stress vector contains NaN or Inf.');
assert(all(isfinite(status1(:))), ...
  'SoftSoilCreepTest:NonFiniteStatus', 'Updated state variables contain NaN or Inf.');
assert(all(isfinite(D(:))), ...
  'SoftSoilCreepTest:NonFiniteTangent', 'Constitutive tangent matrix contains NaN or Inf.');

% 7. Tangent Matrix Symmetry Check
D1 = D(:,:,1);
normD = norm(D1, 'fro');
symmetryError = norm(D1 - D1', 'fro') / max(normD, eps);

assert(symmetryError < 1e-12, ...
  'SoftSoilCreepTest:AsymmetricTangent', ...
  'Constitutive tangent is not symmetric to round-off precision (Error = %.3e).', symmetryError);

% 8. Tangent Matrix Positive Definiteness Check
eigenvalues = eig(D1);
assert(all(eigenvalues > 0), ...
  'SoftSoilCreepTest:NotPositiveDefinite', ...
  'Constitutive tangent matrix is not positive definite.');

% 9. Physical Monotonicity Checks
% (a) Preconsolidation pressure p_c,r must not decrease
assert(status1(1,1) >= status0(1,1) - 1e-8, ...
  'SoftSoilCreepTest:PreconsolidationDecreased', ...
  'Preconsolidation pressure decreased during compressive step.');

% (b) Axial stress magnitude in compression should increase
assert(sigma1(1,3) <= sigma0(1,3), ...
  'SoftSoilCreepTest:StressResponseUnphysical', ...
  'Axial compressive stress did not increase under compressive strain increment.');

% 10. Multi-Gauss Point Array Vectorization Sanity Check
nGauss = 2;
sigmaMulti = repmat(sigma0, nGauss, 1);
depsMulti = repmat(deps, nGauss, 1);
statusMulti0 = mat.initializeStatus(sigmaMulti);

[DMulti, sigmaMulti1, statusMulti1] = mat.getStiffnessMatrix(...
  sigmaMulti, depsMulti, dt, statusMulti0, cellID);

assert(size(DMulti, 3) == nGauss, 'DMulti 3rd dimension must equal nGauss.');
assert(size(sigmaMulti1, 1) == nGauss, 'sigmaMulti1 row count must equal nGauss.');
assert(size(statusMulti1, 1) == nGauss, 'statusMulti1 row count must equal nGauss.');

fprintf('SoftSoilCreep unit test PASSED successfully. (Tangent symmetry error: %.3e)\n', symmetryError);

end
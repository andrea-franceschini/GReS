% Run after initializing GReS and MRST. All quantities use SI units.
% Refine nx to check convergence in space and time together.
nx = 20;
L = 100;
A = 1;
phi = 0.2;
K = 100*milli*darcy;
muW = 5*centi*poise;
muO = 1*centi*poise;
rho = 1000;
dx = L/nx;
PV = phi*A*L;
pIn = 201*barsa;
pOut = 200*barsa;

% Quadratic relative permeabilities, initially oil-filled, no capillarity.
% M is the endpoint mobility ratio (water mobility / oil mobility).
M = muO/muW;
fw = @(s) M*s.^2 ./ (M*s.^2 + (1-s).^2);
dfw = @(s) 2*M*s.*(1-s) ./ (M*s.^2 + (1-s).^2).^2;
mob = @(s) s.^2/muW + (1-s).^2/muO;
sShock = 1/sqrt(1+M);
shockSpeed = fw(sShock)/sShock;  % d(x/L)/d(PVI)
pviBT = 1/shockSpeed;

mesh = structuredMesh(nx,1,1,[0 L],[0 1],[0 1]);
mesh.processGeometry
[x, order] = sort(mesh.cells.center(:,1));
inj = order(1);
prod = order(end);
xFaces = linspace(0,L,nx+1)';

mat = Materials();
mat.addFluid();
mat.addSolid('name',"rock",'cellTags',1);
mat.addPorousRock("rock",'permeability',K,'porosity',phi);

% The end-cell completions provide pure-water inflow and a producer that
% removes both phases at the current upstream fractional flow.
domain = Discretizer('grid',mesh,'materials',mat,'boundaries',Boundaries(mesh));
domain.addPhysicsSolver("TwoPhaseFlow");
flow = domain.getPhysicsSolver("TwoPhaseFlow");
kr = @(s) s.^2;  % GReS calls krO with oil saturation, 1-Sw.
params.wet = struct('viscosity',muW,'compressibility',0,'relPerm',kr,'density',rho);
params.nwet = struct('viscosity',muO,'compressibility',0,'relPerm',kr,'density',rho);
flow.setFluidProperties(params);

wells = PeacemanWells(domain);
radius = min(dx,1)/20;
wells.addWell(inj,'name',"inlet",'dir',"z",'radius',radius, ...
    'sign',1,'value',pIn,'comp',[1 0]);
wells.addWell(prod,'name',"outlet",'dir',"z",'radius',radius, ...
    'sign',-1,'value',pOut);
flow.addWellModel(wells);
WI = [wells.wells.WI];

% Initialize pressure to the all-oil steady pressure field. No pressure
% pin is needed because the well controls already fix the pressure datum.
resistance = (L-dx)/(K*A) + sum(1./WI);
qInitial = (pIn-pOut)/(muO*resistance);
p = zeros(nx,1);
p(order) = pIn - qInitial*muO*(1/WI(1) + (x-x(1))/(K*A));
setState(domain,p,"pressure");
setState(domain,zeros(nx,1),"saturation");

% Bounds follow from min/max total mobility on [0,1].
qMin = (pIn-pOut)/((muW+muO)*resistance);
qMax = (pIn-pOut)/(min(muW,muO)*resistance);
dt = 10*2^floor(log2(0.25*(PV/nx)/(qMax*shockSpeed)));
nSteps = ceil(1.25*pviBT*PV/(qMin*dt));
times = (0:nSteps)*dt;

% Keep every backward-Euler endpoint for the discrete flux integrals.
% DtMin above dt/2 makes a failed step stop instead of silently chopping.
simpar = struct('End',times(end),'DtInit',dt,'DtMax',dt,'DtMin',0.75*dt, ...
    'RelativeTolerance',1e-8,'AbsoluteTolerance',1e-9,'MaxNLIteration',20);
out = OutState('printTimes',times,'solvePrintTimes',1);
solver = NonLinearImplicit('domains',domain, ...
    'simulationparameters',SimulationParameters(simpar),'output',out);
solver.simulationLoop();

%%
results = out.results;
assert(isequal([results.time],times),'Expected all backward-Euler endpoints.');
s = [results.saturation];
p = [results.pressure];
qIn = WI(1)*mob(s(inj,:)).*(pIn-p(inj,:));
qOut = WI(2)*mob(s(prod,:)).*(p(prod,:)-pOut);
waterCut = fw(s(prod,:));
assert(all(qIn > 0 & qOut > 0),'Expected throughflow from inlet to outlet.');

vIn = [0,cumsum(diff(times).*qIn(2:end))];
vOutW = [0,cumsum(diff(times).*qOut(2:end).*waterCut(2:end))];
pvi = vIn/PV;
vol = mesh.cells.volume(:);
waterStored = phi*vol.'*(s-s(:,1));
waterDefect = max(abs(waterStored-vIn+vOutW))/PV;
rateDefect = max(abs(qIn(2:end)-qOut(2:end))./qIn(2:end));
gresLog().log(3,sprintf('Maximum water balance defect / PV: %.3e\n',waterDefect));
gresLog().log(3,sprintf('Maximum relative total-rate imbalance: %.3e\n',rateDefect));
assert(waterDefect < 1e-6 && rateDefect < 1e-6, ...
    'Phase or total-volume balance failed.');
assert(pvi(end) > pviBT,'Simulation ended before analytical breakthrough.');

% The analytical outlet water cut jumps from 0 to fw(sShock). Use half
% that jump as a reproducible threshold for the smeared numerical front.
cutThreshold = 0.5*fw(sShock);
jBT = find(waterCut >= cutThreshold,1);
gresLog().log(3,sprintf('Analytical shock saturation: %.6f\n',sShock));
gresLog().log(3,sprintf('Analytical breakthrough: %.6f PVI\n',pviBT));

if isempty(jBT)
    warning('Numerical water cut did not reach the breakthrough threshold.');
else
    j = jBT-1:jBT;
    pviBTnum = interp1(waterCut(j),pvi(j),cutThreshold);
    assert(abs(pviBTnum-pviBT) < 3e-2)
    gresLog().log(3,sprintf('Numerical breakthrough: %.6f PVI (error %+.3e PVI)\n', ...
        pviBTnum,pviBTnum-pviBT));
end


if gresLog().getVerbosity > 2
figure
tiledlayout(2,2)
for target = [0.35,0.75,1.15]*pviBT
    [~,j] = min(abs(pvi-target));
    ell = L*pvi(j);
    sFace = blProfile(xFaces/ell,sShock,shockSpeed,dfw);

    % Behind the shock x = ell*dfw(Sw). The primitive of Sw with
    % respect to x is x*Sw - ell*fw(Sw).
    % It is zero at the shock and everywhere ahead of it.
    primitive = xFaces.*sFace - ell*fw(sFace);
    sExact = diff(primitive)./diff(xFaces);
    sNum = s(order,j);
    l1Error = sum(vol(order).*abs(sNum-sExact))/sum(vol);
    fprintf('PVI %.6f: mean absolute saturation error %.3e\n',pvi(j),l1Error);

    nexttile
    stairs(xFaces,[sExact;sExact(end)],'k-','LineWidth',1.3)
    hold on
    plot(x,sNum,'b.','MarkerSize',10)
    xlabel('Distance [m]')
    ylabel('Water saturation')
    title(sprintf('%.3f PVI',pvi(j)))
    ylim([0 1])
    legend('Analytical cell average','GReS','Location','best')
end

nexttile
pviPost = linspace(pviBT,pvi(end),200);
sOutlet = blProfile(1./pviPost,sShock,shockSpeed,dfw);
sOutlet(1) = sShock;
plot([0,pviBT,pviPost],[0,0,fw(sOutlet)],'k-','LineWidth',1.3)
hold on
plot(pvi,waterCut,'b-')
xlabel('Pore volumes injected')
ylabel('Outlet water cut')
legend('Analytical','GReS','Location','best')
end

function s = blProfile(xi,sShock,shockSpeed,dfw)
% Invert the rarefaction branch only. The state ahead of the shock is zero.
s = zeros(size(xi));
s(xi == 0) = 1;
fan = find(xi > 0 & xi < shockSpeed);
for j = reshape(fan,1,[])
    speed = min(xi(j),dfw(sShock));
    s(j) = fzero(@(sw) dfw(sw)-speed,[sShock,1]);
end
end
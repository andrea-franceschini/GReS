% inspect convergence of multipliers (affected from instabilities)


% run non conforming model, with multipliers in correct order
patch = 3;
switch patch
    case 1
        Fx = 0; %[kPa]
        Fy = -10;
        nu = 0;
    case {2,3}
        Fx = 10; %[kPa]
        Fy = 0;
        nu = 0.25;
end
% Elastic properties and constitutive tensor
E = 100000;
D = zeros(3);
D([1 5]) = 1-nu;
D([2 4]) = nu;
D(9) = 0.5*(1-2*nu);
D = (E/((1+nu)*(1-2*nu)))*D;
Dmat = D;
stab = 'unstable';
mult = load('mult_fine.dat');
multFine = mult(:,2);
N = 3;
L2 = zeros(N+1,1);
% get multipliers coinciding with correct mesh 
for i = 0:N
    nel = 7*2^i;
    nXs = nel+1;
    rat = 1/2;
    nYs = round(nXs);
    nXm = round(nel*rat+1);
    nYm = round(nXm);
    %%
    [~,ty_conf,~,~] = RunNonConfPatchTest('standard',Fx,Fy,Dmat,patch,stab,nXs,nYs,nXm,nYm);
    skip = round(0.1*nXs);
    laff = [0.5/nel; repmat(1/nel,nXs-2,1); 0.5/nel];
    entries = linspace(1,numel(multFine),nXs);
    multAnal = mult(entries);
    err2 = (multAnal-ty_conf).^2;
    err2 = err2(skip+1:end-skip);
    laff = laff(skip+1:end-skip);
    L2(i+1) = sqrt(sum(err2.*laff));
end
plot(linspace(0,1,nXs),ty_conf)


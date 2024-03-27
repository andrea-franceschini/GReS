function [p,ux,uz,norm_p] = mandelSol(alphan,x,t,c,p0,F,nu,nuU,mu,z)

nx = length(x);
nt = length(t);
nz = length(z);


cellsP = arrayfun(@(alphan) bsxfun(@(x,t) (sin(alphan)/(alphan - ...
    sin(alphan)*cos(alphan)))*(cos(alphan*x/max(x))-cos(alphan))*exp(-alphan^2*c*t/max(x)^2),x',t),alphan,'UniformOutput',false);

cellsUx1 = arrayfun(@(alphan) arrayfun(@(t) ((sin(alphan)*cos(alphan))/(alphan - ...
   sin(alphan)*cos(alphan)))*exp(-alphan^2*c*t/max(x)^2),t),alphan,'UniformOutput',false);

cellsUx2 = arrayfun(@(alphan) bsxfun(@(x,t) ...
    ((cos(alphan)/(alphan - sin(alphan)*cos(alphan)))*sin(alphan*x/max(x))...
    *exp(-alphan^2*c*t/max(x)^2)),x',t),alphan,'UniformOutput',false);

cellsUz = arrayfun(@(alphan) arrayfun(@(t) ((sin(alphan)*cos(alphan))/(alphan - ...
   sin(alphan)*cos(alphan)))*exp(-alphan^2*c*t/max(x)^2),t),alphan,'UniformOutput',false);

%converting cell array into 2S matrices
matP = cell2mat(cellsP);
matUx1 = cell2mat(cellsUx1);
matUx2 = cell2mat(cellsUx2);
matUz = cell2mat(cellsUz);
%reshaping matrices and summation along third direction. 
%Dim 1: position
%Dim 2: time
%Dim 3: Series term

seriesP = sum(reshape(matP,nx,nt,[]),3);
seriesUx1 = sum(reshape(matUx1,1,nt,[]),3);
seriesUx2 = sum(reshape(matUx2,nz,nt,[]),3);
seriesUz = sum(reshape(matUz,1,nt,[]),3);

%calculating solutions
p = 2*p0*seriesP;
ux = x'*((F*nu)/(2*mu*max(x)) - (F*nuU)/(mu*max(x))*seriesUx1)+ F/mu*seriesUx2;
uz = z'*(-(F*(1-nu)/(2*mu*max(x)))+(F*(1-nuU)/(mu*max(x)))*seriesUz);
norm_p = p/p0;

end


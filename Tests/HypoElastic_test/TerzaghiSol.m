function [p,u] = TerzaghiSol(u0,z,t,nm,L,c,pL,biot,gamma,G,nu)

nz = length(z);
m = 0:nm;
nt = length(t);

cellsP = arrayfun(@(m) bsxfun(@(z,t) (1/(2*m+1))*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*sin(((2*m+1)*pi*(L-z))/(2*L)),z',t),m,'UniformOutput',false);
cellsU = arrayfun(@(m) bsxfun(@(z,t) (1/(2*m+1)^2)*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*cos(((2*m+1)*pi*(L-z))/(2*L)),z',t),m,'UniformOutput',false);

matP = cell2mat(cellsP);
matU = cell2mat(cellsU);

seriesP = sum(reshape(matP,nz,nt,[]),3);
seriesU = sum(reshape(matU, nz,nt,[]),3);

p = (4*gamma*pL/pi)*seriesP; %[kPa]
u = repmat(u0',1,nt) + ((1-2*nu)*biot*gamma*pL/2/G/(1-nu))*(repmat(z',1,nt) - (8*L/pi^2)*seriesU); 

end


function [p0,u0] = iniSol(zu,zp,M,pL,Ku,biot,G)
nzp = length(zp);
p0 = zeros(nzp,1);
p0(1:end) = (biot*M*pL)/(Ku+4*G/3);
u0 = arrayfun(@(zu) 1/(Ku+4*G/3)*pL*(zu),zu);
end


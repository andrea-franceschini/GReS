function [p0,u0] = iniSol(z,M,pL,Ku,biot,G)
nz = length(z);
p0 = zeros(nz,1);
p0(1:end) = (biot*M*pL)/(Ku+4*G/3);
u0 = arrayfun(@(z) 1/(Ku+4*G/3)*pL*(z),z);
end


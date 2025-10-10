function [p0,ux0,uz0] = mandelIni(B,F,x,z,nuU,mu)
%calculate initial solution for mandel benchamrk
p0 = (1/(3*max(x)))*B*(1+nuU)*F;
ux0 = F*nuU*x/(2*mu*max(x));
uz0 = -F*(1-nuU)*z/(2*mu*max(x));
end


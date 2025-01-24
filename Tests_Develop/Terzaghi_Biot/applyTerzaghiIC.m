function stateNew = applyTerzaghiIC(stateIn,mat,mesh,bc)
% Apply boundary conditions for Terzaghi problem
p0(1:end) = (biot*M*pL)/(Ku+4*G/3);
u0 = arrayfun(@(zu) 1/(Ku+4*G/3)*pL*(zu),zu);
end


function errNorm = computeInterpError(E,fM,fS,lNod)
% Quadratic error of interpolation for 1D mortar benchmarks
fInt = E * fM;
err2 = (fS - fInt).^2;
errNorm = sqrt(sum(err2.*lNod));
end


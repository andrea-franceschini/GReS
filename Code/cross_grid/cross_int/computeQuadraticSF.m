function N = computeQuadraticSF(pts)
% evaluate quadratic shape functions on coords between -1 and +1
% Return nP x 3 matrix with values of shape functions (to the node)
N(:,1) = -0.5*pts.*(1-pts);
N(:,2) = 1-pts.^2;
N(:,3) = 0.5*pts.*(1+pts);
end

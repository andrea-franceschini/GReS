function N = computeBasis1D(pts,deg)
% coompute 1D basis functions on a set of pointsv in the reference space
switch deg
    case 1 % linear basis functions
        N(:,1) = 0.5 - 0.5*pts;
        N(:,2) = 0.5 + 0.5*pts;
    case 2 % quadratic basis functions
        N(:,1) = -0.5*pts.*(1-pts);
        N(:,2) = 1-pts.^2;
        N(:,3) = 0.5*pts.*(1+pts);
end
end


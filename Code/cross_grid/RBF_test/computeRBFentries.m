function rbf_row = computeRBFentries(dist,type,r)
% given an array of distances, return a row of the RBF matrix
switch type
    case 'wendland'
        rbf_row = pos(1-dist./r).^4.*(1+4*dist./r);
    case 'tps'
        dist(dist == 0) = eps;
        rbf_row = ((dist/r).^2).*log(dist/r);
    case 'gauss'
        rbf_row = exp(-dist.^2/r^2);
    case 'imq'
        rbf_row = (dist.^2+r^2).^(-0.5);
end
end


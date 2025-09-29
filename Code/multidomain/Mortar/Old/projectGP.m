function [xiMaster] = projectGP(coord,topolM,topolS,elM,elS,xi,n,xInt)
% Project a Gauss Point on the slave onto the master side

% Reference: Farah PHD thesis

% elM: master element onto which project
% elS: index of the slave element, origin of the projection
% xi: reference coordinate of the Gauss Points to project
% n: nodal normals
% xInt: real coordinate of integration point 

% OUTPUT
% reference coordinate in the projected slave
% NaN: convergence not reached

% NL param
xiMaster = zeros(length(xi),1);
tol = 1e-8;
itMax = 8;
nodeM = topolM(elM,:);
nodeS = topolS(elS,:);
% evaluate the basis functions on the integration point
for ii = 1:length(xi)
    iter = 0;
    w = 0;
    Nslave = compute1DBasisF(xi(ii));
    Nmaster = compute1DBasisF(xiMaster(ii));
    % evaluate normal at integration point
    rhs1 = [coord(nodeM(1),1) coord(nodeM(2),1); 
        coord(nodeM(1),2) coord(nodeM(2),2)]*Nmaster;
    nS = Nslave'*n([nodeS(1) nodeS(2)],:);
    rhs = rhs1 - w*nS' - xInt(ii);

    % project Gauss Point
    while (norm(rhs,2) > tol) && (iter < itMax)
        % compute Jacobian of the projection
        iter = iter + 1;
        J1 = -0.5*coord(nodeM(1),:) + 0.5*coord(nodeM(2),:);
        J = [J1', -nS'];
        ds = J\(-rhs);
        xiMaster(ii) = xiMaster(ii) + ds(1);
        w = w + ds(2);
        Nmaster = compute1DBasisF(xiMaster(ii));
        % evaluate normal at integration point
        rhs1 = [coord(nodeM(1),1) coord(nodeM(2),1);
            coord(nodeM(1),2) coord(nodeM(2),2)]*Nmaster;
        rhs = rhs1 - w*nS' - xInt(ii,:)';
    end

    if iter >= itMax
        xiMaster(ii) = nan;
        % mark gauss point that did not converge
    end
end


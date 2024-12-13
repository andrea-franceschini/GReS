function [sigma,varargout] = computeStress(mesh, elem, D, u, gauss)
% compute stiffness matrix for a certain domain
sigma = zeros(mesh.nSurfaces,3);
for el = 1:mesh.nSurfaces
   % Get the right material stiffness for each element
   switch mesh.surfaceVTKType(el)
      case 5 % Triangle
         N = getDerBasisF(elem.tri,el);
         B = zeros(3,6);
         B(elem.indB2D(1:12,2)) = N(elem.indB2D(1:12,1));
         u_loc = u(getDoF(mesh.surfaces(el,:)));
         strain = B*u_loc;
         sigma(el,:) = (D*strain)';
      case 9 % Hexahedra % CHECK LATER!
         [N,dJWeighed] = getDerBasisFAndDet(elem.quad,el,1);
         [A,~] = findAreaAndCentroid(elem.quad,el);
         B = zeros(3,8,gauss.nNode);
         B(elem.indB2D(:,2)) = N(elem.indB2D(:,1));
         u_loc = u(getDoF(mesh.surfaces(el,:)));
         strain = pagemtimes(B,u_loc);
         stress_loc = pagemtimes(D,strain);
         stress_loc = stress_loc.*reshape(dJWeighed,1,1,[]);
         sigma(el,:) = sum(stress_loc,3)/A;
   end
end
% nodal projection of stress components

end



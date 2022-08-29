function K = assemblyKGauss(nTotNode, elemMAT, mat, element, varargin)
% Function for global stiffness matrix assembly with no boundary conditions

Kdim = element.dofmax*nTotNode;

%Global stiffness matrix memory allocation
K = sparse(Kdim,Kdim,0);

%Assembly local stiffness matrices for each element
  for el = 1: element.nElem
      k = zeros(element.dofmax*element.nNode);
      
      %Gauss integration points coordinates and weights
      [gCoord,gWeight] = gaussValue(element.nNode,element.dim);
      [~,ncol] = size(gCoord);
      
      % Sum of stiffness matrices calculated in each Gauss points
      % j = n-th Gauss point
      pGauss = 1;
      while pGauss <= ncol
          s = gCoord(1,:);
          t = gCoord(2,:);
          p = gCoord(3,:);
          s = s(pGauss);
          t = t(pGauss);
          p = p(pGauss);
          
          element.getJacobian(el,s,t,p);
          element.getDerivatives(el,s,t,p);
          Bt = (element.B)';
          
          % Get the right material stiffness for each element 
          for i = 1 : length(mat.db.keys)
              if elemMAT(el) == i
                 var = mat.db.keys;
                 matIdentifier = var(i);
                 matIdentifier = char(matIdentifier);
              end
          end
          D = mat.getMaterial(matIdentifier).getStiffnessMatrix(varargin);

          % Calculate stiffness contribution in n-th Gauss point
          for d = 1:element.dim
              W = gWeight(pGauss,d);
              k_pGauss = Bt*D*(element.B)*(det(element.J))*W;
          end

          % Assembly the element local stiffness matrix
          k = k + k_pGauss; 
          pGauss = pGauss+1;
     end

 
       
    dof = length(k)/element.nNode;   
    
    for i_loc = 1:element.nNode
        i_glob = element.elemTopol(el,i_loc);
      
      for j_loc = 1:element.nNode
          j_glob = element.elemTopol(el,j_loc);
          
          ii_global = (i_glob-1)*element.dofmax+1:(i_glob-1)*element.dofmax+3;
          jj_global = (j_glob-1)*element.dofmax+1:(j_glob-1)*element.dofmax+3;
          
          ii_local = (i_loc-1)*dof+1:i_loc*dof;
          jj_local = (j_loc-1)*dof+1:j_loc*dof;

          K(ii_global,jj_global) = K(ii_global,jj_global)+k(ii_local,jj_local);
      end
    end
  end


  for i = 1 : Kdim
       if K(i,i) == 0
           error ('Singular matrix')
       end
  end


end


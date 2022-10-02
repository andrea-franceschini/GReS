function K = assemblyK(nTotNode, elemMAT, mat, element, varargin)
 % Function for global stiffness matrix assembly with no boundary conditions

  % Global stiffness matrix dimension
  Kdim = element.dofmax*nTotNode;

  % Global stiffness matrix memory allocation
  K = sparse(Kdim,Kdim,0);

  for el = 1: element.nElem
    % Get the right material stiffness for each element 
    for i = 1 : length(mat.db.keys)
      if elemMAT(el) == i
         var = mat.db.keys;
         matIdentifier = var(i);
         matIdentifier = char(matIdentifier);
      end
    end
    % Material stiffness matrix
    if nargin > 4
       sz_vec = varargin{1};
       sz = sz_vec(el);
       D = mat.getMaterial(matIdentifier).getStiffnessMatrix(sz);
    else 
       D = mat.getMaterial(matIdentifier).getStiffnessMatrix();
    end
 
    % Assembly local stiffness matrices for each element
    Bt = (element.B(6*el-5:6*el,:))';
    k = element.Vol(el)*Bt*D*element.B(6*el-5:6*el,:); 

    % Number of degrees of freedom for each node
    dof = length(k)/element.nNode;   
    
    % Assembly local stiffness matrix terms inside global stiffness matrix
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

  % Singular matrix check
  for i = 1 : Kdim
     if K(i,i) == 0
         error ('Singular matrix')
     end
  end

end


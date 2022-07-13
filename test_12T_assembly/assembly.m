function [K,F] = assembly(nTotNode, D, element, nodedisp, nodeforce)
% Function for global stiffness matrix assembly with boundary conditions
% and known term

Kdim = element.dofmax*nTotNode;

% 1) Assembly global stiffness matrix

%Global stiffness matrix allocation
K = sparse(Kdim,Kdim,0);

for el = 1: element.nElem
    
    %Assembly local stiffness matrices for each element (è
    %meglio farlo nella classe tetra o nella classe element?)
    Bt = (element.B(6*el-5:6*el,:))';
    k = element.Vol*Bt*D*element.B(6*el-5:6*el,:); 
    
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


% 2) Assembly known term
F = zeros(nodeforce.dofmax*nTotNode,1);
% ciclo switch che distingue il tipo di carico?
for i = 1:Kdim
    for j = 1:length(nodeforce.node_ind)
      if i == nodeforce.node_ind(j)
      F(i) = F(i)+nodeforce.node_value(j);
      end
    end
end


% 3) Assembly boundary conditions (penalty method)
for i = 1:length(nodedisp.node_ind)
    j = nodedisp.node_ind(i);
      F(j) = 0;
      for col = 1:Kdim
          for row = 1:Kdim
              if row == j
                  if row == col
                      K(row,col) = 10^15;
                  end
              end
          end
      end
end
                      

end


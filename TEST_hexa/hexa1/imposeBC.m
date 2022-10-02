function [K,F] = imposeBC(nTotNode, K, dir, neu)
% Function for boundary conditions imposition in global stiffness matrix

  % 1) Known term
  F = zeros(neu.dofmax*nTotNode,1);
  for i = 1:length(K)
    for j = 1:length(neu.bound_dof)
      if i == neu.bound_dof(j)
         F(i) = F(i) + neu.bound_value(j);
      end
    end
  end


  % 2) Imposing boundary conditions (penalty method)
  % Dirichlet value R
  R = max(K, [], 'all');
  for i = 1:length(dir.bound_dof)
    j = dir.bound_dof(i);
    % Known term vector
    F(j) = dir.bound_value(i)*(R*10^10);
    for col = 1:length(K)
       for row = 1:length(K)
          if row == j
              if row == col
                 K(row,col) = R*10^10;
              end
          end
       end
    end
  end
                      
end


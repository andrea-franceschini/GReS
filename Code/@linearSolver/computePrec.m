
% Function for the computation of the preconditioner
function computePrec(obj,A)

   % Check if the problem has multiphysics
   if obj.multiPhysFlag
      % research
      error('not implemented yet');
   else
      % Single physics, single domain, no lagrange multipliers
      if numel(A) == 1 || ~iscell(A)
         % Case of the single Physics preconditioner for one single block
         obj.computeSinglePhPrec(A);

      else
         % Single physics multidomain or lagrange multipliers
         if ~iscell(A)
            error('Passed a non cell matrix to RACP');
         end

         if obj.multidomFlag
            % research
            error('not implemented yet');

         else
            % Treat the multiple domains as if they were one and then use RACP
            if numel(A) ~= 4
               nn = size(A,1);
               A11 = cell2matrix(A(1:nn-1,1:nn-1));
               A12 = cell2matrix(A(1:nn-1,nn));
               A21 = cell2matrix(A(nn,1:nn-1));
               A22 = cell2matrix(A(nn,nn));
   
               clear A;
   
               A = {A11, A12; A21 A22};
            end
            
            % RACP for single physics multi domain 
            obj.computeRACP(A);
         end
      end
   end
end


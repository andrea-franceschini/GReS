
% Function for the computation of the preconditioner
function Compute(obj,A,symMat)

   % Check if it is the case of a single physics single domain
   if numel(A) == 1 || ~iscell(A)
      % Case of the single Physics preconditioner for one single block
      obj.computeSinglePhPrec(A,symMat);
      return
   end

   % Check if the problem has multiphysics
   if obj.multiPhysFlag
%      if numel(A) ~= 4 && iscell(A)
         % research
         error('not implemented yet');
%      else
%         % Single domain Multiphysics
%         obj.computeMCP(A);
%      end
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
            % checks if any of the matrix entries are 0x0 blocks
            [ZeroSpRow,ZeroSpCol] = find(cellfun(@(x) isempty(x), A));

            if ~isempty(ZeroSpRow)
               % get the correct number of rows
               rows = zeros(size(A,1),1);
               for i = 1:size(A,1)
                  rows(i) = size(A{i,1},1);
               end 

               % assign the correct dimension to the matrices
               for i = 1:length(ZeroSpRow)
                  A{ZeroSpRow(i),ZeroSpCol(i)} = sparse(rows(ZeroSpRow(i)),rows(ZeroSpCol(i)));
               end
            end

            nn = size(A,1);
            idxMain = 1:nn-obj.nInt;
            idxSupp = nn+1-obj.nInt:nn;

            A11 = cell2matrix(A(idxMain,idxMain));
            A12 = cell2matrix(A(idxMain,idxSupp));
            A21 = cell2matrix(A(idxSupp,idxMain));
            A22 = cell2matrix(A(idxSupp,idxSupp)); 

            % Fuse the symMat
            symMat1 = zeros(2,2);
            mm = length(idxMain) * length(idxMain);
            sm = length(idxMain) * length(idxSupp);
            ss = length(idxSupp) * length(idxSupp);
            symMat1(1,1) = (sum(symMat(idxMain,idxMain), 'all') == mm);
            symMat1(1,2) = (sum(symMat(idxMain,idxSupp), 'all') == sm);
            symMat1(2,2) = (sum(symMat(idxSupp,idxSupp), 'all') == ss);
            symMat1(2,1) = symMat1(1,2);
            
            clear symMat;
            symMat = symMat1;
            clear symMat1;
            clear A;
   
            A = {A11, A12; A21 A22};
         end

         % RACP for single physics multi domain 
         obj.computeRACP(A,symMat);
      end
   end
end


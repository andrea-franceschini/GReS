% Function to compute the RACP preconditioner for the lagrange multiplier case (single physics multi domain)
function Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   % checks if any of the matrix entries are 0x0 blocks
   [ZeroSpRow,ZeroSpCol] = find(cellfun(@(x) isempty(x), A));
   if ~isempty(ZeroSpRow)
      % get the correct number of rows
      rows = zeros(size(A,1),1);
      for i = 1:size(A,1)
         % Loop over the 
         for j = 1:size(A,1)
            sizz = size(A{i,j},1);
            if sizz ~= 0
               rows(i) = sizz;
               break;
            end
            if j == size(A,1)
               error("no full blocks in this matrix at one column");
            end
         end
      end 
      % assign the correct dimension to the matrices
      for i = 1:length(ZeroSpRow)
         A{ZeroSpRow(i),ZeroSpCol(i)} = sparse(rows(ZeroSpRow(i)),rows(ZeroSpCol(i)));
      end
   end

   nn = size(A,1);
   idxMain = 1:obj.nDom;
   idxSupp = obj.nDom+1:nn;

   % Treat the multiple domains as if they were one and then use RACP
   if numel(A) ~= 4
      
      A11 = cell2matrix(A(idxMain,idxMain));
      A12 = cell2matrix(A(idxMain,idxSupp));
      A21 = cell2matrix(A(idxSupp,idxMain));
      A22 = cell2matrix(A(idxSupp,idxSupp)); 
   
      clear A;
  
      A = {A11, A12; A21 A22};
   end
   
   % Treat the boundary conditions for AMG
   A = obj.treatDirBC(A,obj.PrecSym);
  
   % Compute the augmented matrix
   [A11_aug,inv_D22] = obj.cpt_localAug(A{1,1},A{1,2},A{2,1},A{2,2},obj.PrecSym);
   
   % Compute the test space
   if(obj.phys == 0) % fluids
      TV0 = [];
      for i = 1:obj.nDom - obj.nInt
         TV = ones(size(A{i,1},1),1);
         TV0 = [TV0;TV];
      end
   elseif(obj.phys == 1 || obj.phys == 1.1) % true contact mechanichs physics is 1.1, general poromechanics is 1
      TV0 = [];
      for i = 1:obj.nDom
         TV = mk_rbm_3d(obj.domain(i).grid.coordinates);
         TV0 = [TV0;TV];
      end
   end

   % if obj.DEBUGflag
   %    obj.TV0 = TV0;
   % end
   
   % Compute the amg for block 11
   obj.AMG.Compute(A11_aug,obj.PrecSym,TV0,true);
  
   obj.Apply_L = @(x) obj.ApplyLeft(x,A11_aug,A{1,2},A{2,1},inv_D22);
   obj.Apply_R = @(x) obj.ApplyRight(x);
   
end


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

   % ======================================================================
   % If the user requests that no multidomain treatment is done then
   % compute the racp assuming that there are no multiple domains and or
   % interfaces
   % ======================================================================
   % if obj.multidom == false
   
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
         % if obj.DEBUGflag
         %    obj.TV0 = TV0;
         % end
      end
   
      % Compute the amg for block 11
      obj.AMG.Compute(A11_aug,obj.PrecSym,TV0,true);
   
      obj.Apply_L = @(x) obj.ApplyLeft(x,A11_aug,A{1,2},A{2,1},inv_D22);
      obj.Apply_R = @(x) obj.ApplyRight(x);
   % else
   %    % ===================================================================
   %    % Use the treatment for keeping the multiple domains separate
   %    % ===================================================================
   % 
   %    % If there is more than one interface there is the need to treat it
   %    % carefully
   %    if obj.nInt > 1
   %       Amod = cell(obj.nDom+1,obj.nDom+1);
   %       symMatMod = false(obj.nDom+1,obj.nDom+1);
   % 
   %       % Copy the different physical domains in Amod
   %       Amod(idxMain,idxMain) = A(idxMain,idxMain);
   %       symMatMod(idxMain,idxMain) = symMat(idxMain,idxMain);
   % 
   %       % Loop over the domains to accumulate the interfaces
   %       for i = 1:obj.nDom
   %          % Accumulate the top right interfaces
   %          Amod{i,obj.nDom+1} = cell2matrix(A(i,idxSupp));
   % 
   %          % Accumulate the bottom left interfaces
   %          Amod{obj.nDom+1,i} = cell2matrix(A(idxSupp,i));
   % 
   %          % Accumulate also the symmetry matrix
   %          symMatMod(i,obj.nDom+1) = all(symMat(i,idxSupp),'all');
   %          symMatMod(obj.nDom+1,i) = all(symMat(idxSupp,i),'all');
   %       end
   % 
   %       % Accumulate the lagrange multipliers
   %       Amod{obj.nDom+1,obj.nDom+1} = cell2matrix(A(idxSupp,idxSupp));
   %       symMatMod(obj.nDom+1,obj.nDom+1) = all(symMat(idxSupp,idxSupp),'all');
   % 
   %       clear A symMat;
   % 
   %       A = Amod;
   %       symMat = symMatMod;
   %    end
   % 
   %    % Treat the boundary conditions for AMG
   %    A = obj.treatDirBCmulti(A,symMat);
   % 
   %    % Get the full interfaces matrices
   %    B1 = vertcat(A{1:obj.nDom,end});
   %    B2 = horzcat(A{end,1:obj.nDom});
   %    A11 = cell2matrix(A(idxMain,idxMain));
   % 
   %    % Allocate cell array
   %    A11_aug = cell(obj.nDom,1);
   % 
   %    % Compute the inverse Chat matrix
   %    [inv_D22] = obj.cpt_invD(A11,B1,B2,A{end,end});
   % 
   %    % Treat the domains separately with a block Jacobi approach
   %    for i = 1:obj.nDom
   % 
   %       % Get the submatrix
   %       Ai = [A(i,i), A(i,end); A(end,i), A(end,end)];
   % 
   %       % Compute the addendum
   %       ADD = Ai{1,2}*inv_D22*Ai{2,1}; 
   % 
   %       % Strong Symmetrization if the matrix was seen as symmetric
   %       symm = symMat(i,i) && symMat(i,end) && symMat(end,end);
   %       if symm
   %          ADD = 0.5*(ADD+ADD');
   %       end
   %       A11_aug{i} = Ai{1,1}+ADD;
   % 
   %       % Compute the test space
   %       if(obj.phys == 0) % fluids
   %          TV0 = ones(size(Ai{1,1},1),1);
   %       elseif(obj.phys == 1 || obj.phys == 1.1) % true contact mechanichs physics is 1.1, general poromechanics is 1
   %          TV0 = mk_rbm_3d(obj.domain(i).grid.coordinates);
   %       end
   % 
   %       % Compute the amg for block 11
   %       obj.AMG{i}.Compute(A11_aug{i},symm,TV0,true);
   %       clear TV0;
   %    end
   % 
   %    obj.Apply_L = @(x) obj.ApplyLeft(x,A11_aug,B1,B2,inv_D22);
   %    obj.Apply_R = @(x) obj.ApplyRight(x);
   % end
end


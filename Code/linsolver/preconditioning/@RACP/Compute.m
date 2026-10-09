% Function to compute the RACP preconditioner for the lagrange multiplier case (single physics multi domain)
function A = Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);
   
   % Make the matrix 2x2
   A = obj.condenseDomains(A);

   % Compute and apply Ruiz diagonal scaling if requested
   A = obj.Ruiz.Compute(A);

   % Compute the augmented matrix
   [A11_aug,inv_D22] = obj.cpt_localAug(A{1,1},A{1,2},A{2,1},A{2,2},obj.PrecSym);
   
   % Compute or use provided test space
   if nargin >= 4 && ~isempty(varargin{1})
      TV0 = varargin{1};
   else
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
   end

   TV0 = obj.Ruiz.applyDinv(TV0,1);

   % if obj.DEBUGflag
   %    obj.TV0 = TV0;
   % end
   
   % Compute the inner preconditioner for block 11
   obj.innerPrec.Compute(A11_aug,obj.PrecSym,TV0,true);
  
   obj.Apply_L = @(x) obj.ApplyLeft(x,A11_aug,A{1,2},A{2,1},inv_D22);
   obj.Apply_R = @(x) obj.ApplyRight(x);
   
end


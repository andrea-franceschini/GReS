% Function to compute the Embedded Fractures Preconditioner
function Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   % Treat the boundary conditions in the mechanics block
   A = obj.treatDirBC(A,obj.PrecSym);

   % Compute the test space
   TV0 = [];
   for i = 1:obj.nDom
      TV = mk_rbm_3d(obj.domain(i).grid.coordinates);
      TV0 = [TV0;TV];
   end
   
   % Make A22 slightly more diagonally dominant
   A22 = obj.regularizeA22(A{2,2});
  
   % Compute factorization of the 22 block
   obj.FSAI.Compute(A22,symMat(2,2),true);

   % Compute reverse Schur complement of A11
   invC = @(x) obj.FSAI.ApplyLeft(x);
   obj.S = A{1,1} - A{1,2}*invC(A{2,1});
  
   % Check in case of symmetric indefinite systems for FSAI
   AMGSym = obj.FSAI.PrecSym && obj.PrecSym;

   % Compute the amg for block 11 (mechanics)
   obj.AMG.Compute(obj.S,AMGSym,TV0,true);
   
   % Define handles for application of the preconditioner
   obj.Apply_L = @(x) obj.ApplyLeft(x,A{1,2},A{2,1},invC);
   obj.Apply_R = @(x) obj.ApplyRight(x);

end



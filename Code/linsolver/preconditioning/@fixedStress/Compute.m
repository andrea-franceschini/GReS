% Function to compute the RACP preconditioner for the lagrange multiplier case (single physics multi domain)
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

   % Compute the amg for block 11 (mechanics)
   obj.AMGMech.Compute(A{1,1},symMat(1,1),TV0,true);

   % Compute the approximated Schur complement
   obj.domain.getPhysicsSolver("BiotFullyCoupled").computeRelaxationMatrix();
   S = A{2,2} + (1/obj.problemsolver.dt)*obj.domain.getPhysicsSolver("BiotFullyCoupled").R;

   % Compute the amg for block 22 (fluids)
   obj.AMGFlux.Compute(S,symMat(2,2),ones(size(S,1),1),true);

   obj.Apply_L = @(x) obj.ApplyLeft(x,S,A{1,1},A{1,2});
   obj.Apply_R = @(x) obj.ApplyRight(x);
end


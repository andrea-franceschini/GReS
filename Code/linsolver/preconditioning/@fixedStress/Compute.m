% Function to compute the Fixed Stress preconditioner
function A = Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   % Compute and apply Ruiz diagonal scaling if requested
   A = obj.Ruiz.Compute(A);

   % Compute the test space
   TV0 = [];
   for i = 1:obj.nDom
      TV = mk_rbm_3d(obj.domain(i).grid.coordinates);
      TV0 = [TV0;TV];
   end

   % Scale the test space if needed
   TV0 = obj.Ruiz.applyDinv(TV0,1);

   % Compute the amg for block 11 (mechanics)
   obj.AMGMech.Compute(A{1,1},symMat(1,1),TV0,true);

   % Compute the approximated Schur complement
   obj.domain.getPhysicsSolver("BiotFullyCoupled").computeRelaxationMatrix();
   RR = (1/obj.problemsolver.dt)*obj.domain.getPhysicsSolver("BiotFullyCoupled").R;
   S = A{2,2} + obj.Ruiz.scaleMat(RR,2);

   % Scale the test space if needed
   TVFlux = obj.Ruiz.applyDinv(ones(size(S,1),1),2);

   % Compute the amg for block 22 (fluids)
   obj.AMGFlux.Compute(S,symMat(2,2),TVFlux,true);

   obj.Apply_L = @(x) obj.ApplyLeft(x,S,A{1,1},A{1,2});
   obj.Apply_R = @(x) obj.ApplyRight(x);
end


% Function to compute the Fixed Stress preconditioner
function A = Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   % Determine whether mechanics AMG can be recycled
   recycleMech = obj.isMechLinear && ~isempty(obj.AMGMech.Prec);

   % Compute and apply Ruiz diagonal scaling if requested at the fixedStress level
   if recycleMech
      A = obj.Ruiz.Compute(A, 'freezeBlocks', 1);
   else
      A = obj.Ruiz.Compute(A);

      % Compute the test space (only needed when computing AMGMech)
      TV0 = [];
      for i = 1:obj.nDom
         TV = mk_rbm_3d(obj.domain(i).grid.coordinates);
         TV0 = [TV0;TV];
      end

      % Scale the test space if needed
      TV0 = obj.Ruiz.applyDinv(TV0,1);

      % Store mechanics matrix and compute mechanics AMG
      obj.AMGMech.Compute(A{1,1},symMat(1,1),TV0,true);
   end

   % Time the flow part in case the mechanics is reused
   t_flow_start = tic;

   if strcmpi(obj.flowPrecType,'CPR')
      % Coupled mechanics and two-phase flow (3x3 block system as in White et al. 2019)
      % Block 1: Displacements (u)
      % Block 2: Pressure (p)
      % Block 3: Saturation (s)
      
      % Store coupling block
      Afu = [A{2,1}; A{3,1}];

      % Compute fixed-stress relaxation matrix D^(FS)
      solver = obj.domain(1).getPhysicsSolver("MultiphasePoromechanics");
      solver.computeRelaxationMatrix();
      RR = solver.R;
      
      pState = obj.domain(1).getState("pressure");
      swState = obj.domain(1).getState("saturation");
      cIdx = obj.domain(1).dofm.getFieldCells("pressure");
      tpSolver = solver.getFlowSolver();

      sp = pState(cIdx);
      sw = swState(cIdx);

      rhoW = tpSolver.props.wet.density(sp);
      rhoNW = tpSolver.props.nwet.density(sp);

      dt = obj.problemsolver.dt;

      % Keep diagR sparse; do not convert to full to preserve matrix sparsity
      diagR = diag(RR);
      d11_fs = (1 / dt) * (sw .* rhoW) .* diagR;
      d21_fs = (1 / dt) * ((1 - sw) .* rhoNW) .* diagR;

      % Scale relaxation matrix blocks directly using Ruiz scaling function
      D11_mat = obj.Ruiz.scaleDiag(d11_fs, 2);
      D21_mat = obj.Ruiz.scaleDiag(d21_fs, 3, 2);

      % Form approximate Schur complement flow operator S_ff^(FS) = Aff + D^(FS)
      S_ff = {A{2,2} + D11_mat, A{2,3}; ...
              A{3,2} + D21_mat, A{3,3}};
      
      % Get mechanics dimension for ApplyLeft
      n1 = size(A{1,1}, 1);

      % Compute CPR preconditioner for S_ff
      obj.FlowPrec.Compute(S_ff);
      
      % Bind required blocks directly into Apply_L closure
      obj.Apply_L = @(b) obj.ApplyLeft(b, n1, Afu);
      obj.Apply_R = @(b) obj.ApplyRight(b);

   else
      % Single-phase flow
      % Get the relaxation matrix
      solver = obj.domain.getPhysicsSolver("BiotFullyCoupled");
      solver.computeRelaxationMatrix();
      RR = (1 / obj.problemsolver.dt) * solver.R;
      
      % Compute the approximated Schur complement
      S = A{2,2} + obj.Ruiz.scaleMat(RR,2);

      % Scale the test space if needed
      TVFlux = obj.Ruiz.applyDinv(ones(size(S,1),1),2);

      % Compute the amg for block 22 (fluids)
      obj.FlowPrec.Compute(S,symMat(2,2),TVFlux);

      obj.Apply_L = @(x) obj.ApplyLeft(x,A{1,1},A{1,2});
      obj.Apply_R = @(x) obj.ApplyRight(x);
   end

   obj.tSetupFlow = toc(t_flow_start);
end


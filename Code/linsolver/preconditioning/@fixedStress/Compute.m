% Function to compute the Fixed Stress preconditioner
function A = Compute(obj,A,symMat,varargin)

   if ~iscell(A)
      error("Parent preconditioner not configured to run with multiphys multidomain");
   end

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   % Resolve local block indices
   uBlk = obj.getPhysicsBlocks("displacements");
   pBlk = obj.getPhysicsBlocks("pressure");
   sBlk = obj.getPhysicsBlocks("saturation");
   domId = [obj.blockInfo(uBlk).domainId];

   % Retrieve or compute cached unscaled test space
   if nargin >= 4 && ~isempty(varargin{1})
      TV0 = varargin{1};
   else
      TV0 = obj.getBlockTestSpace(obj.problemsolver);
   end

   % Determine whether mechanics preconditioner can be recycled
   mechLin = obj.isMechLinear();
   recycleMech = mechLin && obj.isMechLinearOld;
   obj.isMechLinearOld = mechLin;

   % Compute and apply Ruiz diagonal scaling if requested at the fixedStress level
   if recycleMech
      A = obj.Ruiz.Compute(A, 'freezeBlocks', uBlk);
   else
      % Scale the matrix and test space if needed
      A = obj.Ruiz.Compute(A);
      TV0 = obj.Ruiz.applyDinv(TV0);

      % Store mechanics matrix and compute mechanics AMG
      TV0_mech = TV0(uBlk);
      obj.MechPrec.Compute(A(uBlk, uBlk), symMat(uBlk, uBlk), TV0_mech, true);
   end

   % Get mechanics part of the matrix and the test spaces
   A_uu = cell2matrix(A(uBlk, uBlk));
   TV0_flow = TV0(pBlk);

   % Time the flow part in case the mechanics is reused
   t_flow_start = tic;

   if strcmpi(obj.flowPrecType,'CPR')
      % Coupled mechanics and two-phase flow (3x3 block system as in White et al. 2019)
      % Store coupling block
      Afu = cell2matrix(A([pBlk, sBlk], uBlk));

      % Compute fixed-stress relaxation matrix D^(FS)
      solver = obj.domain(domId).getPhysicsSolver("MultiphasePoromechanics");
      solver.computeRelaxationMatrix();
      RR = solver.R;
      
      pState = obj.domain(domId).getState("pressure");
      swState = obj.domain(domId).getState("saturation");
      cIdx = obj.domain(domId).dofm.getFieldCells("pressure");
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
      D11_mat = obj.Ruiz.scaleDiag(d11_fs, pBlk);
      D21_mat = obj.Ruiz.scaleDiag(d21_fs, sBlk, pBlk);

      % Form approximate Schur complement flow operator S_ff^(FS) = Aff + D^(FS)
      S_ff = {cell2matrix(A(pBlk,pBlk)) + D11_mat, cell2matrix(A(pBlk,sBlk)); ...
              cell2matrix(A(sBlk,pBlk)) + D21_mat, cell2matrix(A(sBlk,sBlk))};
      
      % Get mechanics dimension for ApplyLeft
      n1 = size(A_uu, 1);

      % Compute CPR preconditioner for S_ff
      obj.FlowPrec.Compute(S_ff, symMat, TV0_flow);
      
      % Bind required blocks directly into Apply_L closure
      obj.Apply_L = @(b) obj.ApplyLeftCPR(b, n1, Afu);
      obj.Apply_R = @(b) obj.ApplyRight(b);

   elseif strcmpi(obj.flowPrecType,'AMG')
      % Single-phase flow
      % Compute relaxation matrix across pressure block(s)
      RR_total = [];
      for b = pBlk
         blk = obj.blockInfo(b);
         if ~contains(string(blk.solverType), "MeshTyingTPFA")
            solver_d = obj.domain(blk.domainId(1)).getPhysicsSolver("BiotFullyCoupled");
            solver_d.computeRelaxationMatrix();
            RR_d = (1 / obj.problemsolver.dt) * solver_d.R;
            RR_total = blkdiag(RR_total, RR_d);
         else
            % MeshTyingTPFA interface pressure: zero relaxation matrix
            RR_total = blkdiag(RR_total, sparse(blk.nDof, blk.nDof));
         end
      end
      
      % Compute the approximated Schur complement
      A_pp = cell2matrix(A(pBlk, pBlk));
      A_up = cell2matrix(A(uBlk, pBlk));
      S = A_pp + obj.Ruiz.scaleMat(RR_total, pBlk);

      % Compute the preconditioner for pressure block
      sym_pp = floor(sum(symMat(pBlk, pBlk), 'all') / (numel(pBlk))^2);
      obj.FlowPrec.Compute(S, sym_pp, TV0_flow);

      obj.Apply_L = @(x) obj.ApplyLeftAMG(x, A_uu, A_up);
      obj.Apply_R = @(x) obj.ApplyRight(x);
   else
      % RACP coupling flow
      % Get the coupling block for pressure and mechanics
      A_up = cell2matrix(A(uBlk, pBlk));

      % Compute relaxation matrix across pressure block(s)
      pBlkvar = intersect(obj.getKindBlocks("variable"),pBlk);
      App = A;
      for b = pBlkvar
         blk = obj.blockInfo(b);
         solver_d = obj.domain(blk.domainId(1)).getPhysicsSolver("BiotFullyCoupled");
         solver_d.computeRelaxationMatrix();
         RR_d = (1 / obj.problemsolver.dt) * solver_d.R;

         % Compute the approximated Schur complement
         App{b,b} = App{b,b} + obj.Ruiz.scaleMat(RR_d, b);
      end

      % Get the correct test space for RACP
      TV0_flow = TV0(pBlkvar);

      % Compute the preconditioner for pressure block
      obj.FlowPrec.Compute(App(pBlk,pBlk), symMat(pBlk,pBlk), TV0_flow);

      obj.Apply_L = @(x) obj.ApplyLeftAMG(x, A_uu, A_up);
      obj.Apply_R = @(x) obj.ApplyRight(x);
   end

   obj.tSetupFlow = toc(t_flow_start);
end
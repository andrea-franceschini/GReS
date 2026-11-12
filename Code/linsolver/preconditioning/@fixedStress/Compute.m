% Function to compute the Fixed Stress preconditioner
function A = Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   % Resolve local block indices
   uBlk = obj.getVarBlocks("displacements");
   pBlk = obj.getVarBlocks("pressure");
   sBlk = obj.getVarBlocks("saturation");
   domId = obj.blockInfo(uBlk(1)).domainId;

   % Determine whether mechanics AMG can be recycled
   recycleMech = obj.isMechLinear && ~isempty(obj.AMGMech.Prec);

   % Compute and apply Ruiz diagonal scaling if requested at the fixedStress level
   if recycleMech
      A = obj.Ruiz.Compute(A, 'freezeBlocks', uBlk);
   else
      % Compute the test space (only needed when computing AMGMech)
      TV0 = obj.AMGMech.buildTestSpace(obj.problemsolver, 1:numel(obj.AMGMech.blockInfo));

      % Scale the matrix and test space if needed
      A = obj.Ruiz.Compute(A);
      TV0 = obj.Ruiz.applyDinv(TV0, uBlk);

      % Store mechanics matrix and compute mechanics AMG
      obj.AMGMech.Compute(A{uBlk,uBlk}, symMat(uBlk,uBlk), TV0, true);
   end

   % Time the flow part in case the mechanics is reused
   t_flow_start = tic;

   if strcmpi(obj.flowPrecType,'CPR')
      % Coupled mechanics and two-phase flow (3x3 block system as in White et al. 2019)
      % Store coupling block
      Afu = [A{pBlk,uBlk}; A{sBlk,uBlk}];

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
      S_ff = {A{pBlk,pBlk} + D11_mat, A{pBlk,sBlk}; ...
              A{sBlk,pBlk} + D21_mat, A{sBlk,sBlk}};
      
      % Get mechanics dimension for ApplyLeft
      n1 = size(A{uBlk,uBlk}, 1);

      % Compute CPR preconditioner for S_ff
      obj.FlowPrec.Compute(S_ff);
      
      % Bind required blocks directly into Apply_L closure
      obj.Apply_L = @(b) obj.ApplyLeft(b, n1, Afu);
      obj.Apply_R = @(b) obj.ApplyRight(b);

   else
      % Single-phase flow
      % Get the relaxation matrix
      solver = obj.domain(domId).getPhysicsSolver("BiotFullyCoupled");
      solver.computeRelaxationMatrix();
      RR = (1 / obj.problemsolver.dt) * solver.R;
      
      % Compute the approximated Schur complement
      S = A{pBlk,pBlk} + obj.Ruiz.scaleMat(RR, pBlk);

      % Scale the test space if needed
      TVFlux = obj.Ruiz.applyDinv(ones(size(S,1),1), pBlk);

      % Compute the amg for pressure block
      obj.FlowPrec.Compute(S, symMat(pBlk,pBlk), TVFlux);

      obj.Apply_L = @(x) obj.ApplyLeft(x, A{uBlk,uBlk}, A{uBlk,pBlk});
      obj.Apply_R = @(x) obj.ApplyRight(x);
   end

   obj.tSetupFlow = toc(t_flow_start);
end


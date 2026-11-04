% Function to compute the Fixed Stress preconditioner
function A = Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   % Determine whether mechanics AMG can be recycled
   recycleMech = obj.isMechLinear && ~isempty(obj.AMGMech.Prec);

   % Compute and apply Ruiz diagonal scaling if requested at the fixedStress level
   if ~isempty(obj.Ruiz) && obj.Ruiz.scalingFlag
      if recycleMech
         A = obj.Ruiz.Compute(A, 'freezeBlocks', 1);
      else
         A = obj.Ruiz.Compute(A);
      end
   end

   % Compute the test space (only needed when computing AMGMech)
   if ~recycleMech
      TV0 = [];
      for i = 1:obj.nDom
         TV = mk_rbm_3d(obj.domain(i).grid.coordinates);
         TV0 = [TV0;TV];
      end

      % Scale the test space if needed
      if ~isempty(obj.Ruiz) && obj.Ruiz.scalingFlag
         TV0 = obj.Ruiz.applyDinv(TV0,1);
      end
   end

   switch upper(obj.flowPrecType)
      case 'CPR'
         % Coupled mechanics and two-phase flow (3x3 block system as in White et al. 2019)
         % Block 1: Displacements (u)
         % Block 2: Pressure (p)
         % Block 3: Saturation (s)
         Auu = A{1,1};
         Aup = A{1,2};
         Aus = A{1,3};
         Apu = A{2,1};
         App = A{2,2};
         Aps = A{2,3};
         Asu = A{3,1};
         Asp = A{3,2};
         Ass = A{3,3};

         % Store mechanics matrix and compute mechanics AMG
         obj.Auu = Auu;
         if ~recycleMech
            obj.AMGMech.Compute(Auu, symMat(1,1), TV0, true);
         end

         % Store coupling blocks
         obj.Auf = [Aup, Aus];
         obj.Afu = [Apu; Asu];

         % Time the flow preconditioner setup (excluding linear mechanics AMG)
         t_flow_start = tic;

         % Compute fixed-stress relaxation matrix D^(FS) (White et al. 2019, Eqs. 24-26)
         dt = obj.problemsolver.dt;
         mSolver = [];
         for d = 1:numel(obj.domain)
            try
               if ismember("MultiphasePoromechanics", obj.domain(d).solverNames)
                  mSolver = obj.domain(d).getPhysicsSolver("MultiphasePoromechanics");
                  break;
               elseif ismember("BiotFullyCoupled", obj.domain(d).solverNames)
                  mSolver = obj.domain(d).getPhysicsSolver("BiotFullyCoupled");
                  break;
               end
            catch
            end
         end

         if ~isempty(mSolver)
            mSolver.computeRelaxationMatrix();
            R = mSolver.R;
         else
            R = sparse(size(App, 1), size(App, 2));
         end

         % Extract saturation and fluid densities for mass accumulation scaling
         try
            pState = obj.domain(1).getState("pressure");
            swState = obj.domain(1).getState("saturation");
            cIdx = obj.domain(1).dofm.getFieldCells("pressure");
            tpSolver = obj.domain(1).getPhysicsSolver("MultiphasePoromechanics").getFlowSolver;
            if isempty(tpSolver) && ~isempty(mSolver)
               try
                  tpSolver = mSolver.getFlowSolver();
               catch
                  tpSolver = [];
               end
            end
            rhoW = full(value(tpSolver.props.wet.density(pState(cIdx))));
            rhoNW = full(value(tpSolver.props.nwet.density(pState(cIdx))));
            sw = full(swState(cIdx));
         catch
            sw = 0.5 * ones(size(App, 1), 1);
            rhoW = ones(size(App, 1), 1);
            rhoNW = ones(size(App, 1), 1);
         end

         diagR = full(diag(R));
         d11_fs = (1 / dt) * (sw .* rhoW) .* diagR;
         d21_fs = (1 / dt) * ((1 - sw) .* rhoNW) .* diagR;

         % Scale relaxation matrix if Ruiz scaling is active on outer blocks
         if ~isempty(obj.Ruiz) && obj.Ruiz.scalingFlag
            D2 = obj.Ruiz.D{2};
            D3 = obj.Ruiz.D{3};
            d11_fs = D2 .* d11_fs .* D2;
            d21_fs = D3 .* d21_fs .* D2;
         end

         nc = size(App, 1);
         D11_mat = sparse(1:nc, 1:nc, d11_fs, nc, nc);
         D21_mat = sparse(1:nc, 1:nc, d21_fs, nc, nc);

         % Form approximate Schur complement flow operator S_ff^(FS) = Aff + D^(FS)
         S_ff = {App + D11_mat, Aps; ...
                 Asp + D21_mat, Ass};
         obj.S_ff = S_ff;

         % Compute CPR preconditioner for S_ff (Ruiz disabled at CPR level!)
         if ~isempty(obj.FlowPrec.Ruiz)
            obj.FlowPrec.Ruiz.scalingFlag = false;
         end
         obj.FlowPrec.Compute(S_ff);
         obj.tSetupFlow = toc(t_flow_start);

         % Set function handles for ApplyLeft
         obj.Apply_L = @(b) obj.ApplyLeft(b);
         obj.Apply_R = @(b) obj.ApplyRight(b);

      otherwise
         % Single-phase flow
         obj.Auu = A{1,1};
         if ~recycleMech
            obj.AMGMech.Compute(A{1,1},symMat(1,1),TV0,true);
         end

         t_flow_start = tic;
         % Compute the approximated Schur complement
         obj.domain.getPhysicsSolver("BiotFullyCoupled").computeRelaxationMatrix();
         RR = (1/obj.problemsolver.dt)*obj.domain.getPhysicsSolver("BiotFullyCoupled").R;
         S = A{2,2} + obj.Ruiz.scaleMat(RR,2);
         obj.S_ff = S;
         obj.Auu = A{1,1};
         obj.Auf = A{1,2};

         % Scale the test space if needed
         TVFlux = obj.Ruiz.applyDinv(ones(size(S,1),1),2);

         % Compute the amg for block 22 (fluids)
         obj.AMGFlux.Compute(S,symMat(2,2),TVFlux,true);
         obj.tSetupFlow = toc(t_flow_start);

         obj.Apply_L = @(x) obj.ApplyLeft(x,S,A{1,1},A{1,2});
         obj.Apply_R = @(x) obj.ApplyRight(x);
   end

end



% Function for the computation of the preconditioner
function computePrec(obj,A)

   % Case with A sparse matrix
   if obj.precOpt == 0
      if iscell(A)
         A = A{1,1};
      end

      if (norm(A-A','f')/norm(A,'f') > 1e-7)
         obj.params.symm = false;
         if DEBUGflag
            fprintf('matrix nonsymmatric %e\n',norm(A-A','f')/norm(A,'f'));
         end
      else
         obj.params.symm = true;
      end

      time_start = tic;
      switch obj.PrecType

         % Compute the AMG preconditioner
         case 'amg'
            
            % Treat Boundary conditions 
            lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);
            A(A==1) = lmax/10;

            % Compute the test space
            if(obj.phys == 0) % fluids
               TV0 = ones(size(A,1),1);
            elseif(obj.phys == 1)
               TV0 = mk_rbm_3d(obj.domain.grid.topology.coordinates);
            end
            
            set_DEBINFO();

            % Actually compute the AMG
            obj.Prec = cpt_aspAMG(obj.params,A,TV0);

            % Define Mfun
            obj.MfunL = @(r) AMG_Vcycle(obj.Prec,A,r);
            obj.MfunR = @(r) r;

         % Compute the FSAI preconditioner
         case 'fsai'
            smootherOp = smoother(A,obj.params.symm,obj.params.smoother);

            % Define Mfun
            [obj.MfunL,obj.MfunR] = defineMfunFSAI(obj,smootherOp);
      end
      T_setup = toc(time_start);
      obj.aTimeComp = obj.aTimeComp + T_setup;
      obj.nComp = obj.nComp + 1;
      
   elseif obj.precOpt == 1

      simple_flag = false;

      % Get the dimensions of the blocks
      n11 = size(A{1,1},1);
      n22 = size(A{1,2},2);

      % If block 22 has dim 0 then resize it
      if size(A{2,2},1) ~= n22
         A{2,2} = sparse(n22,n22);
      end

      % Treat Dirichlet boundary conditions
      A{1,1} = A{1,1}';
      D = sum(spones(A{1,1}));
      ind_dir_dof = find(D==1);
      ind_col_rem = find(sum(spones(A{1,2}))==1);
      [ind_dir_lag,~,~] = find(A{1,2}(:,ind_col_rem));
      ind_dir = union(ind_dir_dof,ind_dir_lag);
      A{1,1}(:,ind_dir) = 0;
      A{1,1} = A{1,1}';
      A{2,1}(:,ind_dir) = 0;
      A{1,2} = A{2,1}';
      fac = max(D);
      D = zeros(n11,1);
      D(ind_dir,1) = fac;
      A{1,1} = A{1,1} + diag(sparse(D));


      % Set RACP Gamma to 1
      gamma = 1.0;

      GLOBAL = false;

      if GLOBAL
         % Compute global augmentation

         lmax_glo = eigs(A{1,1},1,'lm','Display',1,'Tolerance',1.e-5,...
                        'MaxIterations',20,'FailureTreatment','keep');
         fprintf('Eigmax_glo: %15.6e\n',lmax_glo);
      
         ADD_scaled = A{1,2}*A{2,1};
         D = full(diag(ADD_scaled));
         D = D(D>0);
         gmean_ADD_s = geomean(D);
         fprintf('gmean_ADD_s: %15.6e\n',gmean_ADD_s);
      
         ind = find(full(diag(ADD_scaled)) > 0);
         X = A{1,1}(ind,ind);
         fprintf('Gloabl Augmentation factor: %e\n',lmax_glo/gmean_ADD_s);
      
         inv_D22 = gamma*(lmax_glo/gmean_ADD_s)*speye(n22);
      
         A11_aug = A{1,1} + A{1,2}*inv_D22*A{2,1}; A11_aug = 0.5*(A11_aug+A11_aug');
      else

         % Compute local augmentation
         AA_list = {};
         BB_list = {};
         aug = zeros(size(A{2,2},1),1);
         D_11 = full(diag(A{1,1}));
         mean_diag_A = mean(D_11);
         D_22 = full(diag(A{2,2}));
         A21_scaled_T = A{2,1}';
         for icol = 1:n22
            v12 = A{1,2}(:,icol);
            v21 = A21_scaled_T(:,icol);
            [ii_12,~,bb_12] = find(v12);
            [ii_21,~,bb_21] = find(v21);
            if (numel(ii_12)+numel(ii_21) > 0);
               BB = bb_12*bb_21';
               if simple_flag
                  m_a= max(D_11(ii_12));
                  m_b = max(diag(BB));
               else
                  AA = A{1,1}(ii_12,ii_21);
                  AA_list{icol} = AA;
                  BB_list{icol} = BB;
                  m_a = max(eig(full(AA)));
                  m_b = max(eig(full(BB)));
               end
               m_a_sav = m_a;
               if m_a == 0
                  m_a = mean_diag_A;
               end
               alpha = m_a / m_b;
               aug(icol) = 1 / alpha;
            end
         end

         aug_mat = diag(sparse(aug));
         A22_aug = A{2,2} - gamma*aug_mat;

         % Compute augmented 11 block
         inv_D22 = -inv(diag(diag(A22_aug)));
         ADD = A{1,2}*inv_D22*A{2,1}; ADD = 0.5*(ADD+ADD');
         A11_aug = A{1,1}+ADD;
      end

      % Make the computation of the first block preconditioner available 
      obj.precOpt = 0;
      % For now impose the amg
      obj.PrecType = 'amg';

      % Compute the amg for block 11
      computePrec(obj,A11_aug);

      obj.MfunL = @(x) apply_RevAug(obj.Prec,A11_aug,A{1,2},inv_D22,x);
      obj.MfunR = @(x) x;

   else
      error('not implemented yet');
   end
end

































% Helper function for computePrec
function set_DEBINFO()
   global DEBINFO;
   % PARTE GENERALE
   DEBINFO.flag = true;
   DEBINFO.flag = false;
   % PARTE PER LA PROLONGATION
   DEBINFO.prol = [];
   % STAMPARE SI/NO
   DEBINFO.prol.prt_flag = false;
   % UNITA DI STAMPA
   DEBINFO.prol.ofile = 0;
   % STAMPA NUMERO DI ITERAZIONI NEL CALCOLO PROL
   DEBINFO.prol.it_print = false;
   % STAMPA LISTA VICINI NEL CALCOLO PROL
   DEBINFO.prol.neigh_print = false;

   % PARTE PER IL COARSENING
   DEBINFO.coarsen = [];
   % STAMPARE SI/NO
   DEBINFO.coarsen.draw_dist = false;
end







% Helper function for computePrec
% Function to determine how MfunL and MfunR are for each fsai preconditioner
function [MfunL,MfunR] = defineMfunFSAI(obj,smootherOp)
   omega = smootherOp.omega;
   if strcmp(lower(obj.params.smoother.method),'afsai_enh')
      F     = smootherOp.left;
      FT    = smootherOp.right;
      W     = smootherOp.W;
      THETA = smootherOp.THETA;
      MfunL = @(x) omega*(FT*(F*x) + W*(THETA*(W'*x)));
      MfunR = @(x) omega*(FT*(F*x) + W*(THETA*(W'*x)));
   else
      if smootherOp.LS_deg > 0
         MfunL = smootherOp.polyPrec;
         MfunR = smootherOp.polyPrec;
      else
         if strcmp(lower(obj.params.smoother.method),'blk_j')
            MfunL = smootherOp.BLKJ;
            MfunR = smootherOp.BLKJ;
         elseif strcmp(lower(obj.params.smoother.method),'bafsai')
            MfunL = smootherOp.BAFSAI;
            MfunR = smootherOp.BAFSAI;
         elseif strcmp(lower(obj.params.smoother.method),'ddsw')
            MfunL = smootherOp.DDSW1;
            MfunR = smootherOp.DDSW2;
         else
            if (numel(smootherOp.left_out) + numel(smootherOp.right_out)) == 0
               % Simple smoother
               MfunL = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
               MfunR = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
            else
               % Nested smoother
               MfunL = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                                   (smootherOp.left_out*(smootherOp.left*x))));
               MfunR = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                                   (smootherOp.left_out*(smootherOp.left*x))));
            end
         end
      end
   end
end

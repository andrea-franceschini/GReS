% Function to compute the MCP preconditioner for the lagrange multiplier case (multi physics single domain)
function computeMCP(obj,A)

   simple_flag = false;

   % Get the dimensions of the block
   n22 = size(A{1,2},2);

   % Treat Dirichlet boundary conditions
   obj.treatDirBC(A);

   % Compute global augmentation
   lmax_glo = eigs(A{1,1},1,'lm','Display',0,'Tolerance',1.e-5,...
                  'MaxIterations',20,'FailureTreatment','keep');
   
   ADD_scaled = A{1,2}*A{2,1};

   D = full(diag(ADD_scaled));
   D = D(D>0);
   gmean_ADD_s = geomean(D);

   if obj.DEBUGflag
      fprintf('lmax_glo:            %e\n',lmax_glo);
      fprintf('geometric mean BBT:  %e\n',gmean_ADD_s);
      fprintf('Augmentation factor: %e\n',lmax_glo/gmean_ADD_s);
   end
   
   % Compute augmentation block
   gamma = 1.0;
   D22_mat = gamma*(lmax_glo/gmean_ADD_s)*speye(n22);
   
   % Compute new saddle-point system peconditioner
   AA = A{1,1} + A{1,2}*D22_mat*A{2,1}; AA = 0.5*(AA+AA');
   BB = A{1,2};
   CC = A{2,2}; CC = 0.5*(CC+CC');

   if obj.DEBUGflag
      tmp = [AA BB; BB' CC];
      print_SpMat('AUGMENTED_K',AA);
      print_SpMat('AUGMENTED_B',BB);
      print_SpMat('AUGMENTED_C',CC);
      clear tmp;
   end
   
   % For now impose the amg
   obj.PrecType = 'amg';

   % Compute the amg for block 11
   obj.computePrec(AA);

   % Compute preconditioner
   MCP_prec = cpt_MCP2(obj.Prec,AA,BB,CC,obj.maxThreads);

   obj.MfunL = @(x) apply_MCP(MCP_prec,AA,BB,CC,x);
   obj.MfunR = @(x) x;
end



% Function to compute the true MCP preconditioner
function MCP_prec = cpt_MCP2(AMG_prec,AA,BB,CC,nthread)

nstep = 30;
step_size = 1;
epsilon = 1.e-3;
G11 = afsai_cpp(AA,nthread,nstep,step_size,epsilon);

HH = G11*BB;
SS = CC - HH'*HH;

G22 = afsai_cpp(-SS,nthread,nstep,step_size,epsilon);

MCP_prec.AMG_prec = AMG_prec;
MCP_prec.G22 = G22;

end

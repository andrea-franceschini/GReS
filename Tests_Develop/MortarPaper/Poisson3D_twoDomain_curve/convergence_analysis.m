clear
close all
addpath(genpath(pwd))

if exist("logConvAnalysis","file") == 2
  delete("logConvAnalysis")
end
diary logConvAnalysis

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

%% Poisson problem with single domain in 3D. Testing new poisson module

% analytical solution
anal = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradx = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2 - 2*X(1) + pi*cos(pi*X(1)));
grady = @(X) -pi*sin(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradz = @(X) -pi*cos(pi*X(2)).*sin(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
f = @(X) cos(pi*X(2)).*cos(pi*X(3)).*h(X(1));


%% model commons
simParam = SimulationParameters('simParam.dat');
% base structure to write xml file
strDomain = readstruct('Domains/domain1block.xml');
interfFile = fullfile('Domains','interfaces_1.xml');
strInterf = readstruct(interfFile);
out_str = []; % output structure
%% INPUT
% base domain size
N_0_l = 2;
N_0_r = 3;

% number of refinement
nref = 6;
[h,L2,H1] = deal(zeros(nref,1));

% study parameters
elem_type = "hexa";                 % hexa,hexa27
integration_type = ["SegmentBased","RBF","ElementBased"];    % SegmentBased (7 gp),ElementBased,RBF

% add one more evenutally

N_l = [4 8 16 32 40 48];

N_r = [6 12 24 48 60 72];

%% convergence loop

for i_t = integration_type
  if strcmp(i_t,'SegmentBased')
    nG = 7;
  else
    nG = [3 6 16];
  end

  nInt = 0;

  if strcmp(i_t,'RBF')
    nInt = 5;
  end
  for ngp = nG
    fprintf('_____________________________________________________________________\n')
    fprintf('Running convergence analysis with %s integration - %i GP \n',i_t,ngp)
    fprintf('_____________________________________________________________________\n')
    for n_i = nInt
      if strcmp(i_t,'RBF')
        fprintf('Using %i interpolation points %i \n',n_i)
      end
      % refinement loop
      for i = 1:nref

        N_i_l = N_l(i);
        N_i_r = N_r(i);
        fprintf('_____________________________________________________________________\n')
        fprintf('Running mesh refinement %i - %s %i GP \n',i,i_t,ngp);
        fprintf('_____________________________________________________________________\n')

        % run script to get refined mesh
        fname = strcat('domain_',elem_type,'_',num2str(i));
        command = "python Mesh/scripts/domain_curve.py "  + fname...
          + " " + num2str(N_i_l) + " " + num2str(N_i_r) + " " + elem_type;
        %system(command);

        %return

%         continue

        meshFile = fullfile('Mesh','meshes',fname+".vtk");
        mesh = Mesh();
        mesh.importMesh(meshFile);

        % set up bc file
        nodes = unique(mesh.surfaces(ismember(mesh.surfaceTag,[2 4]),:));
        c = mesh.coordinates(nodes,:);
        vals = arrayfun(@(i) anal(c(i,:)),1:numel(nodes));
        vals = reshape(vals,[],1);
        writeBCfiles('BCs/bc','NodeBC','Dir','Poisson','manufactured_bc',0,0,nodes,vals);

        clear mesh

        % write mesh to domain file
        strDomain.Domain.Name = "Cube_"+integration_type+"_"+num2str(i);
        strDomain.Domain.Geometry = fullfile('Mesh','meshes',fname+".vtk");
        domainFile = fullfile('Domains','domain1block.xml');
        writestruct(strDomain,domainFile);


        % write interface to file
        strInterf.Interface(1).Quadrature.typeAttribute = i_t;
        strInterf.Interface(1).Quadrature.nGPAttribute = ngp;
        %strInterf.Interface(1).Print.nameAttribute = "interf_"+i_t+"_"+num2str(i);
        if strcmp(i_t,'RBF')
          strInterf.Interface(1).Quadrature.nIntAttribute = n_i;
        end
        writestruct(strInterf,interfFile);

        % processing Poisson problem
        domains = buildModelStruct_new(domainFile,simParam);
        domains.Discretizer.getSolver('Poisson').setAnalSolution(anal,f,gradx,grady,gradz);

        [interfaces,domains] = Mortar.buildInterfaceStruct(interfFile,domains);
        % set up analytical solution

        solver = MultidomainFCSolver(simParam,domains,interfaces);
        solver.NonLinearLoop();
        %solver.finalizeOutput();

        %runPoisson;

        pois = getSolver(domains.Discretizer,'Poisson');
        [L2(i),H1(i)] = pois.computeError_v2();
        h(i) = 1/N_i_r;
        fprintf('Max absolute error is: %1.6e \n',max(abs(pois.state.data.err)));

      end % end refinement

      % compute convergence order
      L2ord = log(L2(1:end-1)./L2(2:end))./log(h(1:end-1)./h(2:end));
      H1ord = log(H1(1:end-1)./H1(2:end))./log(h(1:end-1)./h(2:end));

      % store result in structure
      switch elem_type
        case 'hexa'
          outDir = "Output/OUT_HEXA_2";
        case 'hexa27'
          outDir = "Output/OUT_HEXA27_2";
      end
      out_str = struct('int_type',i_t,'nGP',ngp,...
        'L2norm',L2,'H1norm',H1);
      fname = i_t+"_"+num2str(ngp);
      if strcmp(i_t,'RBF')
        fname = fname + "_" + num2str(n_i);
        out_str.nInt = n_i;
      end
      save(fullfile(outDir,strcat(fname,'.mat')),'-struct',"out_str");
    end % end interpolation loop points
  end % end gp loop
end % end integration type loop

diary off
% store output structure

%% plotting convergence profiles
% figure(1)
% loglog(h,L2,'-ro')
% hold on
% loglog(h,H1,'-b^')
% xlabel('h')
% ylabel('error_norm')
% legend('L2','H1')





% 3d interpolation of an analytical function from one interface to the
% other


clear
close all
% simply define 3D meshes made of 1 line of 3D elements on either sides of
% the interface, compute the mortar operator ( possibly dual), and evaluate
% the error as done in the past

nref = 6;
elem_type = 'hexa27';
N0 = 2;

nInt = 4;


% analytical function to test
f = @(x,y) sin(4*x).*cos(4*y);
L2 = zeros(nref,1);

rbftype = 'gauss';
%
% for i = 1:nref
%   % write mesh to file
%   % set refinement
%   Nb = 1.5*(2*2^(i-1));
%   Nt = (2*2^(i-1));
% 
%   fnameBot = strcat('bot_',num2str(i),'_',elem_type);
%   command = "python Mesh/scripts/domain_bot.py "  + fnameBot...
%     + " " + num2str(Nb) + " " + elem_type;
%   system(command);
% 
%   fnameTop = strcat('top_',num2str(i),'_',elem_type);
%   command = "python Mesh/scripts/domain_top.py "  + fnameTop...
%     + " " + num2str(Nt) + " " + elem_type;
%   system(command);
% 
% end

for i_t = ["RBF"]
  for i = 1:nref

    fnameBot = strcat('bot_',num2str(i),'_',elem_type);
    fnameTop = strcat('top_',num2str(i),'_',elem_type);


    % define models
    strDomain = readstruct('Domains/domain2block.xml');
    strDomain.Domain(1).Name = "top"+num2str(i);
    strDomain.Domain(1).Geometry = fullfile('Mesh','meshes',fnameTop+".vtk");
    strDomain.Domain(2).Name = "bot"+num2str(i);
    strDomain.Domain(2).Geometry = fullfile('Mesh','meshes',fnameBot+".vtk");
    domainFile = fullfile('Domains','domain2block.xml');
    writestruct(strDomain,domainFile);

    if strcmp(i_t,'SegmentBased')
      nG = 7;
    else
      nG = 4;
    end

    % define interface
    strInterf = readstruct(fullfile('Domains','interface.xml'));
    strInterf.Interface(1).Quadrature.typeAttribute = i_t;
    strInterf.Interface(1).Quadrature.nGPAttribute = nG;
    strInterf.Interface(1).Quadrature.nIntAttribute = nInt;
    strInterf.Interface(1).Quadrature.rbfAttribute = rbftype;

    % create instance of Mortar class
    domains = buildModelStruct_new(domainFile,[]);

    mortar = MeshGlueDual(1,strInterf.Interface(1),domains);

    

    mortar.computeMortarMatrices();
    E = mortar.computeMortarOperator;

    % extract only local E
    nm = mortar.mesh.local2glob{1};
    ns = mortar.mesh.local2glob{2};
    E = E(ns,nm);

    % compute analytical function on the master side
    xM = mortar.mesh.msh(1).coordinates(:,1);
    yM = mortar.mesh.msh(1).coordinates(:,2);
    xS = mortar.mesh.msh(2).coordinates(:,1);
    yS = mortar.mesh.msh(2).coordinates(:,2);
    fMaster = f(xM,yM);
    fSlave = E*fMaster;
    L2(i) = computeError(fSlave,f,mortar);
  end

  fname = i_t;
  if strcmp(i_t,'RBF')
    fname = fname + "_" + rbftype + "_" + num2str(nInt);
  end
  switch elem_type
    case 'hexa'
      outname = fullfile('OUT','HEXA',fname+".mat");
    case 'hexa27'
      outname = fullfile('OUT','HEXA27',fname+".mat");
  end
  save(outname,"L2");
end


function L2err = computeError(fAppr,fEx,mortar)
msh = mortar.mesh.msh(2);
L2err = 0;
for el = 1:msh.nSurfaces
  elem = mortar.getElem(2,el);
  N = getBasisFinGPoints(elem);
  [dJW] = getDerBasisFAndDet(elem,el);
  c_gp = getGPointsLocation(elem,el);
  dof = msh.surfaces(el,:);
  f_gp = N*fAppr(dof);
  err = f_gp - fEx(c_gp(:,1),c_gp(:,2));
  err_2 = err.^2;        % squared value of error
  L2errLoc = sum(err_2.*reshape(dJW,[],1));
  L2err = L2err + L2errLoc;
end
L2err = sqrt(L2err);
end



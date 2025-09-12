classdef ElementBasedQuadrature < handle
  % Implement utilities to perform element based integration given a pair
  % of master and slave element in contact
  % Can be equipped with the RBF scheme to cheaply sort out elements that
  % are not really connected

  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    msh       % instance of InterfaceMesh class
    mortar
    nGP
  end

  properties (Access = private)
    % temporary properties for element-based sort out process
    idSlave = 0       % current id of slave element
    tempGPloc         % full list of gp location
    tempdJw           % full list of jacobian determinant
    tempNs            % full list of slave basis functions
    tempNmult         % full list of multiplier basis functions
    tempNbubble       % full list of bubble basis function
    gptoProj          % mark gp still to be projected for current m/s couple
    flagProj          % projected gp for current master/slave couple
  end

  methods
    function obj = ElementBasedQuadrature(mortar,nG)
      obj.mortar = mortar;
      obj.msh = obj.mortar.mesh.msh;
      obj.nGP = nG;
    end

    function [Ns,Nm,Nmult,NbubSlave,NbubMaster] = getMortarBasisFunctions(obj,is,im)

      elemSlave = obj.mortar.getElem(2,is);
      elemMaster = obj.mortar.getElem(1,im);

      % new slave element
      if obj.idSlave ~= is
        if ~isempty(obj.gptoProj)
          if any(obj.gptoProj)
            warning('% i GP not projected for element %i',...
              sum(obj.gptoProj),obj.idSlave)
          end
        end
        obj.idSlave = is;
        obj.tempdJw = getDerBasisFAndDet(elemSlave,is);
        obj.tempNs = getBasisFinGPoints(elemSlave);
        obj.tempNmult = computeMultiplierBasisF(obj.mortar,is,obj.tempNs);
        obj.tempGPloc = elemSlave.GaussPts.coord;
        obj.gptoProj = true(size(obj.tempNs,1),1);
        if nargout > 3
          obj.tempNbubble = getBubbleBasisFinGPoints(elemSlave);
        end
      else         % old slave element, sort out gp
        if ~any(obj.gptoProj)
          % gp already projected on all previous elements
          [Ns,Nm,Nmult,NbubSlave,NbubMaster] = deal([]);
          return
        end
% 
%         obj.tempNs = obj.tempNs(~obj.suppFlag,:);
%         obj.tempNmult = obj.tempNmult(~obj.suppFlag,:);
%         if nargout > 3
%           obj.tempNbubble = obj.tempNbubble(~obj.suppFlag,:);
%         end
%         obj.tempGPloc = obj.tempGPloc(~obj.suppFlag,:);
%         obj.tempdJw = obj.tempdJw(~obj.suppFlag);
      end

      % reset list of projected gp
      obj.flagProj = false(size(obj.tempNs,1),1);
      % get master basis and gp in slave support
      xiProj = projectGP(obj,is,im);
      % remove projected gp from list
      obj.gptoProj(obj.flagProj,:) = false;
      Nm = elemMaster.computeBasisF(xiProj);

      if nargout > 4
        NbubMaster = elemMaster.computeBubbleBasisF(xiProj);  
        NbubMaster = NbubMaster(obj.flagProj,:);
      end

      if all(obj.gptoProj)
        % no detected intersection
        [Ns,Nm,Nmult,NbubSlave,NbubMaster] = deal([]);
        return
      end

      %Nm = Nm(obj.suppFlag,:);

      % return basis function in active gp
      Ns = obj.tempNs(obj.flagProj,:);
      Nmult = obj.tempNmult(obj.flagProj,:);
      if nargout > 3
        NbubSlave = obj.tempNbubble(obj.flagProj,:);
      end
    end

%     function [xiM, gpProjected] = par_projectGP(obj, is, im)
%       elM = getElem(obj.mortar,1,im);
%       elS = getElem(obj.mortar,2,is);
%       nodeS = obj.msh(2).surfaces(is,:);
%       coordS = obj.msh(2).coordinates(nodeS,:);
%       X = elS.Nref(obj.gptoProj,:) * coordS;
%       xiS = obj.tempGPloc(obj.gptoProj,:);
%       ngp = size(xiS,1);
%       xiM_all = zeros(ngp, 2);       % results
%       fl = false(ngp,1);             % flags
%       tol = 1e-9;
%       itMax = 8;
%       k = find(obj.gptoProj);        % global indices of GP to project
% 
%       nodeM = obj.msh(1).surfaces(im,:);
%       coordM = obj.msh(1).coordinates(nodeM,:);
%       normals = obj.mortar.mesh.avgNodNormal{2}(nodeS,:);
% 
%       % Create local copies of elements (safe for parallel)
%       elMcopy = elM;
%       elScopy = elS;
% 
%       parfor i = 1:ngp
%         xiMi = elScopy.centroid;
%         xiSi = xiS(i,:);
%         Xi = X(i,:);
%         Ns = elScopy.computeBasisF(xiSi);
%         ng = Ns * normals;
%         ng = ng / norm(ng);
% 
%         iter = 0;
%         w = 0;
%         Nm = elMcopy.computeBasisF(xiMi);
%         rhs = (Nm * coordM - w * ng - Xi)';
% 
%         while norm(rhs,2) > tol && iter < itMax
%           iter = iter + 1;
%           dN = elMcopy.computeDerBasisF(xiMi);
%           J1 = dN * coordM;
%           J = [J1' -ng'];
%           ds = J \ (-rhs);
%           xiMi = xiMi + ds(1:2)';
%           w = w + ds(3);
%           Nm = elMcopy.computeBasisF(xiMi);
%           rhs = (Nm * coordM - w * ng - Xi)';
%         end
% 
%         if FEM.checkInRange(elMcopy, xiMi)
%           xiM_all(i,:) = xiMi;
%           fl(i) = true;
%         end
%       end
% 
%       xiM = xiM_all(fl,:);
%       gpProjected = fl;
% 
%       % Safe to update obj.flagProj *after* the loop
%       obj.flagProj(k(fl)) = true;
%     end


    function [xiM,gpProjected] = projectGP(obj,is,im)
      % xi: reference coordinates of the gauss point
      % get nodal normal
      elM = getElem(obj.mortar,1,im);
      elS = getElem(obj.mortar,2,is);
      nodeS = obj.msh(2).surfaces(is,:);
      coordS = obj.msh(2).coordinates(nodeS,:);
      X = elS.Nref(obj.gptoProj,:)*coordS;                    % real position of gauss pts
      xiS = obj.tempGPloc(obj.gptoProj,:);
      ngp = size(xiS,1);
      xiM = zeros(ngp,2);
      itMax = 8;
      tol = 1e-14;
      gpProjected = false(size(obj.gptoProj,1),1);
      fl = false(ngp,1);
      k = find(obj.gptoProj); % index of gp to project
      for i = 1:ngp
        Ns = elS.computeBasisF(xiS(i,:));
        ng = Ns*obj.mortar.mesh.avgNodNormal{2}(nodeS,:); % slave normal at GP
        ng = ng/norm(ng);
        xiM(i,:) = elS.centroid;                          % initial guess
        iter = 0;
        w = 0;
        nodeM = obj.msh(1).surfaces(im,:);
        coordM = obj.msh(1).coordinates(nodeM,:);
        Nm = elM.computeBasisF(xiM(i,:));
        rhs = (Nm*coordM - w*ng - X(i,:))';
        %
        while (norm(rhs,2) > tol) && (iter < itMax)
          iter = iter+1;
          dN = elM.computeDerBasisF(xiM(i,:));
          J1 = dN*coordM;
          J = [J1' -ng'];
          ds = J\(-rhs);
          xiM(i,:) = xiM(i,:) + (ds(1:2))';
          w = w + ds(3);
          Nm =  elM.computeBasisF(xiM(i,:));
          rhs = (Nm*coordM - w*ng - X(i,:))';
        end
        % check if gp lies in master reference space and update gp flag 
        if FEM.checkInRange(elM,xiM(i,:))
          % turn off gp that fall in master 
          % turn on flag for projected gp
          obj.flagProj(k(i)) = true;
          fl(i) = true;
          % store result in projected gp
        end
      end
      xiM = xiM(fl,:);
    end


    function mat = integrate(obj,func,varargin)
      % element based integration
      % check input
      assert(nargin(func)==numel(varargin),['Number of specified input (%i)' ...
        'not matching the integrand input (%i)'],numel(varargin),nargin(func));
      size3 = cellfun(@(x) size(x, 3), varargin);
      if ~all(size3 == size3(1))
        error('All inputs must have the same size along dimension 3.');
      end
      dJWeighed = obj.tempdJw(obj.flagProj);
      mat = func(varargin{:});
      mat = mat.*reshape(dJWeighed,1,1,[]);
      mat = sum(mat,3);
    end
    
  end
end
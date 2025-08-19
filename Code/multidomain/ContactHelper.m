classdef ContactHelper < handle
  
  % General helper for contact algorithms  
  % Local rotation matrix computation
  % activeSet update

  properties
    activeSet
    tol      % tolerances for active set check
  end

  properties (Access = private)
    rotationMat             % store global rotation matrices for each element
    normals                 % element normal/slave normals   
    node_normals
    elem
  end

 
  methods
    function obj = ContactHelper(mortar,varargin)
      obj.elem = mortar.elements(2);
      computeNormals(obj,mortar);
      computeRotationMatrix(obj);
      defineTolerances(obj);
      initializeActiveSet(obj,mortar)
    end

    %     function activeSet = updateActiveSet(obj,activeSet)
    %
    %
    %     end

    function initializeActiveSet(obj,mortar)
      % 1 - stick
      % 2 - new slip
      % 3 - slip
      % 4 - open
      nMult = mortar.mesh.msh(2).nSurfaces;
      obj.activeSet = ones(nMult,1);      % all stick at the beginning
    end

    function R = getRotationMatrix(obj,in)
      % get rotation matrix of entity i
      if isscalar(in)
        R = obj.rotationMat(in,:);
        R = reshape(R,3,3);
      else
        % input is a 3D/4D array of normals in each GP
        sz = size(in);
        R = zeros(3,3,sz(3),sz(4));
        for i = 1:sz(4)
          for j = 1:sz(3)
            R(:,:,j,i) = obj.computeRot(in(:,:,j,i));
          end
        end
      end
    end

    function n = computeNormalinGP(obj,i)

      % return a 3D matrix 3x1xnG with the normal evaluated in each gp
      el = getSurfElementByID(obj.elem,i);
      nG = el.GaussPts.nNode;
      n = zeros(3,nG);
      for i = 1:nG
        xg = el.GaussPts.coord(i,:);
        n(:,i) = el.computeNormal(i,xg);
      end
      n = reshape(3,1,nG);
    end

    function computeNormals(obj,mortar)
      % compute normal of each element in the contact interface
      nS = mortar.mesh.msh(2).nSurfaces;
      switch mortar.multiplierType
        case 'P0'
          obj.normals = zeros(nS,3);
          for elId = 1:nS
            elType = getSurfElementByID(obj.elem,elId);
            obj.normals(elId,:) = computeNormal(elType,elId);
          end
        otherwise
          error('Contact solver is already available with P0 multipliers')
      end
      obj.node_normals = mortar.mesh.avgNodNormal{2};
    end

    function n = getNodalNormal(obj,elemId)
      nodeId = obj.elem.mesh.surfaces(elemId,:);
      n = obj.node_normals(nodeId,:);
      n = reshape(n',[],1);
    end

    function n = getNormal(obj,elemId)
      n = obj.normals(elemId,:);
      n = reshape(n,3,1);
    end


    function computeRotationMatrix(obj)
      % compute rotation matrix of each element in the contact interface
      % this maps traction dofs to global reference frame

      nS = size(obj.normals,1);
      obj.rotationMat = zeros(nS,9);

      for i = 1:nS
        n = obj.normals(i,:);
        R = obj.computeRot(n);
        obj.rotationMat(i,:) = R(:);
      end
    end

    function out = isStick(obj,elID)
      out = obj.activeSet(elID)==1;
    end

    function out = isOpen(obj,elID)
      out = obj.activeSet(elID)==4;
    end

    function out = isSlip(obj,elID)
      out = obj.activeSet(elID)==3;
    end

    function out = isNewSlip(obj,elID)
      out = obj.activeSet(elID)==2;
    end

    function setSlip(obj,elID)
      obj.activeSet(elID) = 3;
    end

    function setNewSlip(obj,elID)
      obj.activeSet(elID) = 2;
    end

    function setStick(obj,elID)
      obj.activeSet(elID) = 1;
    end

    function setOpen(obj,elID)
      obj.activeSet(elID) = 4;
    end

    function defineTolerances(obj)
      % for now, we use default values
      obj.tol.sliding = 1e-4;
      obj.tol.normalGap = 1e-5;
      obj.tol.normalTrac = 1e-3;
      obj.tol.slidingCheck = 3e-2;
      obj.tol.areaTol = 1e-2;
    end
  end

  methods (Static)
    function R = computeRot(n)
      n = reshape(n,1,[]);
      % candidate orthogonal basis guess
      m1 = [n(3),0,-n(1)];
      m2 = [0,n(3),-n(2)];
      norm_m1 = norm(m1);
      norm_m2 = norm(m2);

      if norm_m1 + eps > norm_m2
        m2 = cross(n,m1);
        m1 = m1/norm(m1);
        m2 = m2/norm(m2);
      else
        m1 = cross(n,m2);
        m1 = -m1/norm(m1);
        m2 = m2/norm(m2);
      end

      R = [n',m1',m2'];

      assert(abs(det(R)-1.0)<10*eps,'Rotation matrix non unit determinant')
    end
  end

end


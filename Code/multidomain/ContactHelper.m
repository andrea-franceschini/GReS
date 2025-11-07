classdef ContactHelper < handle
  
  % General helper for contact algorithms  
  % Local rotation matrix computation
  % activeSet update

  properties
    activeSet
    tol      % tolerances for active set check
    forceStickBoundary = false
    resetActiveSet = false;
  end

  properties (Access = private)
    rotationMat             % store global rotation matrices for each element
    normals                 % element normal/slave normals   
    node_normals
    elem
    multType

  end

 
  methods
    function obj = ContactHelper(mortar,varargin)
      obj.elem = mortar.elements(2);
      obj.node_normals = mortar.mesh.avgNodNormal{2};
      obj.multType = mortar.multiplierType;

      if strcmp(obj.multType,'P0')
        computeNormals(obj,mortar);
        computeRotationMatrix(obj);
      end

      if strcmp(class(mortar),"ContactMortar")
        defineTolerances(obj);
        initializeActiveSet(obj,mortar)
      end
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
      obj.activeSet.curr = ones(nMult,1);      % all stick at the beginning
      obj.activeSet.prev = obj.activeSet.curr;
      obj.activeSet.stateChange = zeros(nMult,1); % count how many times an element has changed state during iteration
    end

    function R = getRotationMatrix(obj,elemId)
      % get rotation matrix of entity i
      if strcmp(obj.multType,"P0")
        R = obj.rotationMat(elemId,:);
        R = reshape(R,3,3);
      else
        n = getNodalNormal(obj,elemId);
        sz = size(n,1);
        R = zeros(sz(1));
        for i = 1:3:sz(1)
          R(i:i+2,i:i+2) = obj.computeRot(n(i:i+2));
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
      end
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
      out = obj.activeSet.curr(elID)==1;
    end

    function out = isOpen(obj,elID)
      out = obj.activeSet.curr(elID)==4;
    end

    function out = isSlip(obj,elID)
      out = obj.activeSet.curr(elID)==3;
    end

    function out = isNewSlip(obj,elID)
      out = obj.activeSet.curr(elID)==2;
    end

    function setSlip(obj,elID)
      obj.activeSet.curr(elID) = 3;
    end

    function setNewSlip(obj,elID)
      obj.activeSet.curr(elID) = 2;
    end

    function setStick(obj,elID)
      if obj.activeSet.stateChange(elID) < obj.tol.maxStateChange
        % avoid repeated oscillations
        obj.activeSet.curr(elID) = 1;
      end
    end

    function setOpen(obj,elID)
      obj.activeSet.curr(elID) = 4;
    end

    function defineTolerances(obj)
      % for now, we use default values
      obj.tol.sliding = 1e-4;
      obj.tol.normalGap = 1e-6;
      obj.tol.normalTrac = 1e-3;
      obj.tol.slidingCheck = 1e-2;
      obj.tol.minLimitTraction = 1e-4;  % below this value, the limit traction is set to 0
      obj.tol.areaTol = 1e-2;
      obj.tol.maxStateChange = 4;
    end
  end

  methods (Static)
    %     function R = computeRot(n)
    %       n = reshape(n,1,[]);
    %       % candidate orthogonal basis guess
    %       m1 = [n(3),0,-n(1)];
    %       m2 = [0,n(3),-n(2)];
    %       norm_m1 = norm(m1);
    %       norm_m2 = norm(m2);
    %
    %       if norm_m1 + 1e2*eps > norm_m2
    %         m2 = cross(n,m1);
    %         m1 = m1/norm(m1);
    %         m2 = m2/norm(m2);
    %       else
    %         m1 = cross(n,m2);
    %         m1 = -m1/norm(m1);
    %         m2 = m2/norm(m2);
    %       end
    %
    %       R = [n',m1',m2'];
    %
    %       assert(abs(det(R)-1.0)<1e-8,'Rotation matrix non unit determinant')
    %     end


    function R = computeRot(n)

      n = reshape(n,1,[]);
      n = n / norm(n);   % normalize input normal

      % Pick a vector not parallel to n (to start Gram–Schmidt)
      if abs(n(1)) < 0.9
        tmp = [1,0,0];
      else
        tmp = [0,1,0];
      end

      % First tangent: orthogonalize tmp against n
      m1 = tmp - dot(tmp,n)*n;
      m1 = m1 / norm(m1);

      % Second tangent: orthogonal to both
      m2 = cross(n,m1);
      m2 = m2 / norm(m2);

      % Assemble rotation matrix
      R = [n', m1', m2'];

      % Check orientation: enforce det=+1 (right-handed)
      if det(R) < 0
        m1 = -m1;
        R = [n', m1', m2'];
      end

      assert(abs(det(R)-1.0) < 1e-12, ...
        'Rotation matrix not orthogonal to machine precision');
    end
  end

end


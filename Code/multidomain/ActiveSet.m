classdef ActiveSet < handle
  
  % General helper for contact algorithms  
  % Local rotation matrix computation
  % activeSet update

  properties

    state   % active set flag 
    tol     
  end


  methods

    function obj = ActiveSet(interf,varargin)

        defineTolerances(obj);
        initializeActiveSet(obj,interf)
      end

    function initializeActiveSet(interface,N)
      % 1 - stick
      % 2 - new slip
      % 3 - slip
      % 4 - open
      % the interface MUST have an activeSet property
      try
        % all stick at the beginning
      interface.activeSet.curr = repmat(ContactMode.stick,N,1);      
      interface.activeSet.prev = obj.activeSet.curr;
      % count how many times an element has changed state during iteration
      interface.activeSet.stateChange = zeros(N,1); 
      catch
        error("activeSet not found in the interface solver")
      end
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


end


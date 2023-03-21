classdef TransvElastic < handle
  % ELASTIC TRANSVERSE ISOTROPIC material class

  properties (Access = private)
    % Elastic modulus in vertical direction (z), i.e., perpendicular to the
    % symmetry plane
    Ev
    % Poisson's ratio in vertical direction (relates vertical and
    % horizontal strain) (z)
    vv
    % Ratio of horizontal to vertical elastic moduli
    beta     % beta = Eh/Ev
    % Ratio of horizontal to vertical Poisson ratios
    gamma    % gamma = vh/vv
    % Ratio of horizontal to vertical shear moduli
    eta      % eta = Gh/Gv
    %
    %
    % Derived properties
    Gh
    Eh
    vh
    Gv
    % M factor (sigmax/sigmaz)
    M
    % Vertical compressibility Cm
    cM
  end

  methods (Access = public)
    % Class constructor method
    function obj = TransvElastic(inputString)
      % Calling the function to set object properties 
      obj.setMaterialParameters(inputString);
    end

    % Material stiffness matrix calculation using the object properties
    function D = getStiffnessMatrix(obj)
      % Stiffness matrix
      D = zeros(6);
      a = 1/obj.Eh;
      b = -obj.vh/obj.Eh;
      c = -obj.vv/obj.Eh;
      d = 1/obj.Ev;
      D(1:3,1:3) = inv([a b c; b a c; c c d]);
      D(22) = obj.Gh;
      D([29 36]) = obj.Gv;
    end
    
    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end
    
    % Get vertical compressibility
    function cM = getRockCompressibility(obj)
      cM = obj.cM;
    end
  end

  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function setMaterialParameters(obj, block)
      % Preliminary check on the number of rows in each material block
      % and the number of parameters
      nEntry = size(block,1);
      if nEntry ~= 2
        error('Material %s\n Wrong number of input rows', block(2));
      end
      strParams = strsplit(block(2));
      nEntry = size(strParams,2);
      if nEntry ~= 5
        error('Material %s\n Wrong number of input parameters',block(1));
      end
      params = str2double(strParams);
      % Assign object properties
      obj.Ev = params(1);
      obj.vv = params(2);
      obj.beta = params(3);
      obj.gamma = params(4);
      obj.eta = params(5);
      %
      % Compute derived material properties
      obj.Eh = obj.beta * obj.Ev;
      obj.vh = obj.gamma * obj.vv;
      obj.Gh = obj.Eh/(2*(1+obj.vh));
      obj.Gv = obj.Gh/obj.eta;
      %
      % Check if the material parameters fulfill the thermodynamic
      % consistency requirements (see Janna et al, 2012 -
      % doi:10.1016/j.ijrmms.2012.01.015)
      assert(obj.Eh > 0,'The horizontal elastic modulus Eh MUST be POSITIVE');
      assert(obj.Ev > 0,'The vertical elastic modulus Ev MUST be POSITIVE');
      assert(obj.Gv > 0,'The vertical shear modulus Gv MUST be POSITIVE');
      assert((1-(obj.vh)^2) > 0,'1-vh^2 MUST be POSITIVE');
      assert((1-obj.vh-2/(obj.beta)*(obj.vv)^2) > 0,'1-vh-2*Ev/Eh*vv^2 MUST be POSITIVE');
      %
      % Compute the M factor
      obj.M = obj.vv/(1-obj.vv*obj.gamma);
      %
      % Compute vertical compressibility
      obj.cM = 1/obj.Ev*(1-2*(obj.vv)^2/((1-obj.vh)*obj.beta));
    end
  end
end
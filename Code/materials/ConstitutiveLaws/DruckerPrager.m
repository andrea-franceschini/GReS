classdef DruckerPrager < handle
  % ELASTIC ISOTROPIC material class

  properties (Access = public)
    % Elastic modulus
    E
    % Poisson's ratio
    nu
    % M factor (sigmax/sigmaz)
    M
    % Vertical compressibility Cm
    cM
    % Dilatancy angle
    psi
    % Friction angle
    phi
    % Cohestion
    co
    % Hardening parameter
    h
    % Initial hardening variable
    k
    % Drucker Prager coefficients
    alpha
    beta  
    epsilon
    varepsilon
    % flag for tabular input prameters
    isTabular
  end

  methods (Access = public)
     % Class constructor method
     function obj = DruckerPrager(varargin)
        % Calling the function to set the object properties
           obj.readMaterialParameters(varargin{:});
     end

    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,2);
    end
    %
    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status, cellID)
      nptGauss = size(sigmaIn,1);
      D = getElasticTensor(obj,cellID);
      DAll = repmat(D,[1, 1, nptGauss]);
      % trial elastic stress
      sigmaOut0 = sigmaIn + epsilon*D;
      % Update stress and constitutive tensore with Return mapping algorithm
      [sigmaOut, DAll] = returnMapping(obj, sigmaOut0, nptGauss, DAll);
    end
    %

    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end

    % Get vertical compressibility
    function cM = getRockCompressibility(obj)
      cM = obj.cM;
    end
    function [sigma, D] = returnMapping(obj, sigmaIn, nptGauss, D)
      sigma = sigmaIn;
      for i = 1:nptGauss
        q = sqrt(0.5*(((sigmaIn(i, 1)-sigmaIn(i, 2))^2+(sigmaIn(i, 1)- ...
            sigmaIn(i, 3))^2+(sigmaIn(i, 2)-sigmaIn(i, 3))^2))+ ...
            3*(sigmaIn(i, 4)^2+sigmaIn(i, 5)^2+sigmaIn(i, 6)^2));
        qtr = q;
        p = sum(sigmaIn(i, 1:3))/3;
        I = [1;1;1;0;0;0];
        n = ((1.5/q)*(sigmaIn(i, 1:6)'- p*I));
        f = q/sqrt(3) + obj.alpha*p-obj.epsilon*(obj.co+obj.k);
        G = obj.E/(2*(1+obj.nu));
        K = obj.E/(3*(1-2*obj.nu));
        lambdac = max(0, f/(G+obj.alpha*obj.beta*K+(obj.epsilon^2)*obj.h));
        if f > 0
           p = p-lambdac*K*obj.beta;
           q = q-lambdac*sqrt(3)*G;
           obj.varepsilon = obj.varepsilon+lambdac*obj.epsilon;
           if q < 0 %apex return
              sigma(i, 1:6) = (p*I)';
              D(:,:,i) = K*((1-(obj.alpha*obj.beta*K)/(obj.alpha*obj.beta*K+obj.epsilon^2*obj.h)))*(I*I');
           else
              sigma(i, 1:6) = (p*I + 2/3*q*n)';
              var1 = lambdac*(2*sqrt(3)*G^2)/(qtr);
              var2 = diag([1 1 1 0.5 0.5 0.5])-(1/3)*(I*I')-(2/3)*(n*n');
              var3 = (2*G)/(sqrt(3))*n+K*obj.beta*I;
              var4 = (2*G)/(sqrt(3))*n+obj.alpha*K*I;
              var5 = G+obj.alpha*obj.beta*K+obj.epsilon^2*obj.h;
              D(:,:,i) = D(:,:,i) - var1*var2-var3*((var4)/(var5))';
           end
        else
           continue
        end
      end
    end


    function D = getElasticTensor(obj,cID)
       D = zeros(6);
       if obj.isTabular
          pois = obj.nu(cID);
          D([1 8 15]) = 1-pois;
          D([2 3 7 9 13 14]) = pois;
          D([22 29 36]) = (1-2*pois)/2;
          D = obj.E(cID)/((1+pois)*(1-2*pois))*D;
       else
          D([1 8 15]) = 1-obj.nu;
          D([2 3 7 9 13 14]) = obj.nu;
          D([22 29 36]) = (1-2*obj.nu)/2;
          D = obj.E/((1+obj.nu)*(1-2*obj.nu))*D;
       end
    end
  end

  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function readMaterialParameters(obj, varargin)

      default = struct('youngModulus',[],...
                       'poissonRatio',[],...
                       'dilatancy',[],...
                       'cohesion',[],...
                       'hardeningParameter',[],...
                       'hardeningVariable',0.0);

      params = readInput(default,varargin{:});

      obj.E = params.youngModulus;
      obj.nu = params.poissonRatio;
      obj.psi = params.dilatancy;
      obj.phi = params.frictionAngle;
      obj.co = params.cohesion;
      obj.h = params.hardeningParameter;
      obj.varepsilon = params.hardeningVariable;

      % Compute the M factor
      obj.M = obj.nu/(1-obj.nu);
      % Compute vertical compressibility
      obj.cM = (1+obj.nu)*(1-2*obj.nu)/(obj.E*(1-obj.nu));

  
      obj.alpha = (3*tan(deg2rad(obj.phi)))/(sqrt(9+12*tan(deg2rad(obj.phi))^2)); 
      obj.beta = (3*tan(deg2rad(obj.psi)))/(sqrt(9+12*tan(deg2rad(obj.psi))^2));
      obj.epsilon = 3/(sqrt(9+12*tan(deg2rad(obj.phi))^2));
      obj.k = obj.h*obj.varepsilon;

    end



    function readTabMaterialParameters(obj,fID,fileName,mesh)
       % read parameters in tabular format 
       youngModFile = readToken(fID,fileName);
       poissonRatioFile = readToken(fID,fileName);
       obj.E = setTabularParams(youngModFile,mesh);
       obj.nu = setTabularParams(poissonRatioFile,mesh);
       % Compute the M factor
       obj.M = obj.nu./(1-obj.nu);
       % Compute vertical compressibility
       obj.cM = (1+obj.nu).*(1-2*obj.nu)./(obj.E.*(1-obj.nu));
       psiFile = readToken(fID,fileName);
       phiFile = readToken(fID,fileName);
       cohesionFile = readToken(fID,fileName);
       hardParamFile = readToken(fID,fileName);
       hardVarFile = readToken(fID,fileName);
       obj.psi = setTabularParams(psiFile,mesh);
       obj.phi = setTabularParams(phiFile,mesh);
       obj.co = setTabularParams(cohesionFile,mesh);
       obj.h = setTabularParams(hardParamFile,mesh);
       obj.varepsilon = setTabularParams(hardVarFile,mesh);
       obj.alpha = (3*tan(deg2rad(obj.phi)))./(sqrt(9+12*tan(deg2rad(obj.phi)).^2));
       obj.beta = (3*tan(deg2rad(obj.psi)))./(sqrt(9+12*tan(deg2rad(obj.psi)).^2));
       obj.epsilon = 3./(sqrt(9+12*tan(deg2rad(obj.phi)).^2));
       obj.k = obj.h.*obj.varepsilon;
    end
  end
end

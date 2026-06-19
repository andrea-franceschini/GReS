classdef DruckerPrager < handle
  % DRUCKER-PRAGER elastoplastic material class

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
    % Cohesion
    co
    % Hardening parameter
    h
    % Drucker Prager coefficients
    alpha
    beta  
    epsilon
    % flag for tabular input parameters
    isTabular
  end

  methods (Access = public)
     % Class constructor method
     function obj = DruckerPrager(varargin)
        obj.readMaterialParameters(varargin{:});
     end

    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      % Returns 6 columns to remain compatible with Poromechanics state allocation
      status = zeros(nptGauss,6); 
    end

    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilonIn, dt, status, cellID)
      nptGauss = size(sigmaIn,1);
      D = getElasticTensor(obj,cellID);
      DAll = repmat(D,[1, 1, nptGauss]);
      
      % Validated matrix layout for incoming engineering strain increment matrix
      sigmaOut0 = sigmaIn + epsilonIn * D;
      [sigmaOut, DAll, status] = returnMapping(obj, sigmaOut0, nptGauss, DAll, status, cellID);
    end

    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end

    % Get vertical compressibility
    function cM = getRockCompressibility(obj)
      cM = obj.cM;
    end

    function [sigma, D, status] = returnMapping(obj, sigmaIn, nptGauss, D, status, cellID)
      sigma = sigmaIn;
      I = [1;1;1;0;0;0];
      P_dev = diag([1, 1, 1, 0.5, 0.5, 0.5]) - (1/3) * (I * I');
      
      if obj.isTabular
         local_E       = obj.E(cellID);
         local_nu      = obj.nu(cellID);
         local_alpha   = obj.alpha(cellID);
         local_beta    = obj.beta(cellID);
         local_epsilon = obj.epsilon(cellID);
         local_co      = obj.co(cellID);
         local_h       = obj.h(cellID);
      else
         local_E       = obj.E;
         local_nu      = obj.nu;
         local_alpha   = obj.alpha;
         local_beta    = obj.beta;
         local_epsilon = obj.epsilon;
         local_co      = obj.co;
         local_h       = obj.h;
      end

      G = local_E / (2 * (1 + local_nu));
      K = local_E / (3 * (1 - 2 * local_nu));
      
      for i = 1:nptGauss
        p_tr = sum(sigmaIn(i, 1:3)) / 3;
        s_vec = sigmaIn(i, 1:6)' - p_tr * I;
        
        % Closed parentheses syntax fix enclosing shear stress terms
        q_tr = sqrt(0.5 * ((sigmaIn(i,1) - sigmaIn(i,2))^2 + ...
                           (sigmaIn(i,1) - sigmaIn(i,3))^2 + ...
                           (sigmaIn(i,2) - sigmaIn(i,3))^2) + ...
                    3 * (sigmaIn(i,4)^2 + sigmaIn(i,5)^2 + sigmaIn(i,6)^2));
        
        local_varepsilon = status(i,1);
        f = q_tr / sqrt(3) + local_alpha * p_tr - local_epsilon * (local_co + local_h * local_varepsilon);
        
        if f > 0
           denom = G + local_alpha * local_beta * K + (local_epsilon^2) * local_h;
           lambdac = f / denom;
           
           p_new = p_tr - lambdac * K * local_beta;
           q_new = q_tr - lambdac * sqrt(3) * G;
           
           if q_new < 0 || q_tr < 1e-12
              denom_apex = local_alpha * K * local_beta + (local_epsilon^2) * local_h;
              f_trial_apex = local_alpha * p_tr - local_epsilon * (local_co + local_h * local_varepsilon);
              lambda_apex = f_trial_apex / denom_apex;
              
              local_varepsilon = local_varepsilon + lambda_apex * local_epsilon;
              p_apex = p_tr - lambda_apex * K * local_beta;
              
              sigma(i, 1:6) = (p_apex * I)';
              D(:,:,i) = ((local_epsilon^2 * local_h * K) / denom_apex) * (I * I');
           else
              sigma(i, 1:6) = (p_new * I + (q_new / q_tr) * s_vec)';
              local_varepsilon = local_varepsilon + lambdac * local_epsilon;
              
              % Closed-form exact consistent tangent matrix execution
              theta_vec = local_alpha * K * I + (sqrt(3) * G / q_tr) * s_vec;
              psi_vec = local_beta * K * I + (sqrt(3) * G / q_tr) * s_vec;
              
              D(:,:,i) = K * (I * I') + 2 * G * (q_new / q_tr) * P_dev + ...
                         (3 * G / (q_tr^2)) * (1 - q_new / q_tr) * (s_vec * s_vec') - ...
                         (1 / denom) * (psi_vec * theta_vec');
           end
           status(i, 1) = local_varepsilon;
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
    function readMaterialParameters(obj, varargin)
      default = struct('youngModulus',[],...
                       'poissonRatio',[],...
                       'frictionAngle',[],...
                       'dilatancy',[],...
                       'cohesion',[],...
                       'hardeningParameter',[],...
                       'hardeningVariable',0.0);

      params = readInput(default,varargin{:});

      obj.E = params.youngModulus;
      obj.nu = params.poissonRatio;
      obj.phi = params.frictionAngle;
      obj.psi = params.dilatancy;
      obj.co = params.cohesion;
      obj.h = params.hardeningParameter;

      if numel(obj.E) > 1
         obj.isTabular = true;
      else
         obj.isTabular = false;
      end

      obj.M = obj.nu ./ (1 - obj.nu);
      obj.cM = (1 + obj.nu) .* (1 - 2 * obj.nu) ./ (obj.E .* (1 - obj.nu));

      obj.alpha = (3 * tan(deg2rad(obj.phi))) ./ (sqrt(9 + 12 * tan(deg2rad(obj.phi)).^2)); 
      obj.beta = (3 * tan(deg2rad(obj.psi))) ./ (sqrt(9 + 12 * tan(deg2rad(obj.psi)).^2));
      obj.epsilon = 3 ./ (sqrt(9 + 12 * tan(deg2rad(obj.phi)).^2));
    end

    function readTabMaterialParameters(obj,fID,fileName,mesh)
       youngModFile = readToken(fID,fileName);
       poissonRatioFile = readToken(fID,fileName);
       obj.E = setTabularParams(youngModFile,mesh);
       obj.nu = setTabularParams(poissonRatioFile,mesh);
       
       obj.M = obj.nu./(1-obj.nu);
       obj.cM = (1+obj.nu).*(1-2*obj.nu)./(obj.E.*(1-obj.nu));
       psiFile = readToken(fID,fileName);
       phiFile = readToken(fID,fileName);
       cohesionFile = readToken(fID,fileName);
       hardParamFile = readToken(fID,fileName);
       obj.psi = setTabularParams(psiFile,mesh);
       obj.phi = setTabularParams(phiFile,mesh);
       obj.co = setTabularParams(cohesionFile,mesh);
       obj.h = setTabularParams(hardParamFile,mesh);
       
       obj.alpha = (3*tan(deg2rad(obj.phi)))./(sqrt(9+12*tan(deg2rad(obj.phi)).^2));
       obj.beta = (3*tan(deg2rad(obj.psi)))./(sqrt(9+12*tan(deg2rad(obj.psi)).^2));
       obj.epsilon = 3./(sqrt(9+12*tan(deg2rad(obj.phi)).^2));
    end
  end
end
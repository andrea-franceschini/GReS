classdef DKPlastic < handle
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
    psi
    phi
    co
    h
    alpha
    beta
    epsilon
  end

  methods (Access = public)
    % Class constructor method
    function obj = DKPlastic(fID, matFileName)
      % Calling the function to set the object properties 
      obj.readMaterialParameters(fID, matFileName);
    end
    
    function [status] = initializeStatus(obj, sigma)
      nptGauss = size(sigma,1);
      status = zeros(nptGauss,2);
    end
    %
    % Material stiffness matrix calculation using the object properties
    function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status)
      nptGauss = size(sigmaIn,1);
%       sigmaOut = zeros(nptGauss,6);    % MODIFICA SN
      % Stiffness matrix
%       DAll = zeros(6,6,nptGauss);      % MODIFICA SN
      D = zeros(6);
      D([1 8 15]) = 1-obj.nu;
      D([2 3 7 9 13 14]) = obj.nu;
      D([22 29 36]) = (1-2*obj.nu)/2;
      D = obj.E/((1+obj.nu)*(1-2*obj.nu))*D;
      
      % MODIFICA SN
%       for i = 1 : nptGauss
%         sigmaOut(i,:) = sigmaIn(i,:) + epsilon(i,:)*D;
%         DAll(:,:,i) = D;
%       end
      sigmaOut0 = sigmaIn + epsilon*D;
      sigmaOut = lambdacorr(obj, sigmaOut0, nptGauss);
      %fprintf('%+.5E \n',sigmaOut-sigmaOut0);
      DAll = repmat(D,[1, 1, nptGauss]);
    end
    %
    % Material stiffness matrix calculation using the object properties
%     function D = getStiffnessMatrix(obj)
%       % Stiffness matrix
%       D = zeros(6);
%       D([1 8 15]) = 1-obj.nu;
%       D([2 3 7 9 13 14]) = obj.nu;
%       D([22 29 36]) = (1-2*obj.nu)/2;
%       D = obj.E/((1+obj.nu)*(1-2*obj.nu))*D;
%     end
    
    % Method that returns the M factor
    function m = getMFactor(obj)
      m = obj.M;
    end
    
    % Get vertical compressibility
    function cM = getRockCompressibility(obj)
      cM = obj.cM;
    end

    function sigma = lambdacorr(obj, sigmaIn, nptGauss)
              for i = 1:nptGauss
                  q = sqrt(0.5*(((sigmaIn(i, 1)-sigmaIn(i, 2))^2+(sigmaIn(i, 1)- ...
                      sigmaIn(i, 3))^2+(sigmaIn(i, 2)-sigmaIn(i, 3))^2))+ ...
                      sigmaIn(i, 4)^2+sigmaIn(i, 5)^2+sigmaIn(i, 6)^2);
                  p = sum(sigmaIn(i, 1:3)); 
                  f = q/sqrt(3) + obj.alpha*p-obj.epsilon*(obj.co);
                  G = obj.E/(2*(1+obj.nu));
                  K = obj.E/(3*(1-2*obj.nu)); 
                  lambdac = max(0, f/(G+obj.alpha*obj.beta*K+(obj.epsilon^2)*obj.h));
                  fprintf('initial lambdac %+.5E \n', lambdac);
                  sigma = sigmaIn;   
                  while f > 0
                      sigma = sigmaIn;
                      p = p-lambdac*K*obj.beta;
                      q = q-lambdac*sqrt(3)*G;
                      if q < 0 %q <0 e non ha maniera di correggersi
                          f = obj.alpha*p - obj.epsilon*(obj.co);
                          lambdac = max(0, f/(G+obj.alpha*obj.beta*K+(obj.epsilon^2)*obj.h));fprintf('q < 0 and f is %+.5E | lambda is %+.5E | p is %+.5E \n ', f, lambdac, p);
                          sigma(i, 1:6) = (p-lambdac*K*obj.beta).*ones(1,6);    
                      else
                          f = q/sqrt(3) + obj.alpha*p-obj.epsilon*obj.co; 
                          n = (1.5/q).*(sigma(i, 1:6)- p/3.*[1 1 1 0 0 0]);
                          lambdac = max(0, f/(G+obj.alpha*obj.beta*K+(obj.epsilon^2)*obj.h));fprintf('q < 0 and f is %+.5E | lambda is %+.5E | p is %+.5E \n ', f, lambdac, p);
                          sigma(i, 1:6) = sigma(i, 1:6) - lambdac.*((2*G/sqrt(3)).*n+K*obj.beta);                          
                      end
                  end
               end
        end
  end

  methods (Access = private)
    % Assigning material parameters (check also the Materials class)
    % to object properties
    function readMaterialParameters(obj, fID, matFileName)
      tmpVec = readDataInLine(fID, matFileName, 6);
      %
      obj.E = tmpVec(1);
      obj.nu = tmpVec(2);
      %
      % Compute the M factor
      obj.M = obj.nu/(1-obj.nu);
      %
      % Compute vertical compressibility
      obj.cM = (1+obj.nu)*(1-2*obj.nu)/(obj.E*(1-obj.nu));
      obj.psi = tmpVec(3);
      obj.phi = tmpVec(4);
      obj.co = tmpVec(5); %coesione
      obj.h = tmpVec(6); %hardening parameter
      obj.alpha = (3*tan(obj.phi))/(sqrt(9+12*tan(obj.phi)^2)); 
      obj.beta = (3*tan(obj.psi))/(sqrt(9+12*tan(obj.psi)^2));
      obj.epsilon = 3/(sqrt(9+12*tan(obj.phi)^2));
    end
  end
end
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
      DAll = repmat(D,[1, 1, nptGauss]);
      sigmaOut0 = sigmaIn + epsilon*D;
      [sigmaOut, DAll] = lambdacorr(obj, sigmaOut0, nptGauss, DAll);
      
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

    function [sigma, D] = lambdacorr(obj, sigmaIn, nptGauss, D)
              sigma = sigmaIn;
              for i = 1:nptGauss
                  q = sqrt(0.5*(((sigmaIn(i, 1)-sigmaIn(i, 2))^2+(sigmaIn(i, 1)- ...
                      sigmaIn(i, 3))^2+(sigmaIn(i, 2)-sigmaIn(i, 3))^2))+ ...
                      3*(sigmaIn(i, 4)^2+sigmaIn(i, 5)^2+sigmaIn(i, 6)^2));
                  p = sum(sigmaIn(i, 1:3));
                  f = q/sqrt(3) + obj.alpha*p-obj.epsilon*(obj.co);
                  G = obj.E/(2*(1+obj.nu));
                  K = obj.E/(3*(1-2*obj.nu)); 
                  lambdac = max(0, f/(G+obj.alpha*obj.beta*K+(obj.epsilon^2)*obj.h));
                  I = [1 1 1 0 0 0];
                  if f > 0 
                      p = p-lambdac*K*obj.beta;
                      q = q-lambdac*sqrt(3)*G;
                      if q < 0 %apex return
                          q = 0;
                          f = obj.alpha*p - obj.epsilon*(obj.co);
                          lambdac = max(0, f/(obj.alpha*obj.beta*K+(obj.epsilon^2)*obj.h));
                          sigma(i, 1:6) = p.*I;
                          p = p-lambdac*K*obj.beta; %%
                          f = obj.alpha*p - obj.epsilon*(obj.co); %%
                          %sigma(i, 1:6) = p.*I;
                          D(:,:,i) = K*((1-(obj.alpha*obj.beta*K)/(obj.alpha*obj.beta*K+obj.epsilon^2*obj.h))).*(I'*I);
                      else
                          f = q/sqrt(3) + obj.alpha*p-obj.epsilon*obj.co; 
                          n = ((1.5/q).*(sigma(i, 1:6)- p/3.*I))';
                          lambdac = max(0, f/(G+obj.alpha*obj.beta*K+(obj.epsilon^2)*obj.h));
                          sigma(i, 1:6) = sigma(i, 1:6) - lambdac.*((2*G/sqrt(3)).*n'+K*obj.beta.*I);
                          var1 = lambdac*(2*sqrt(3)*G^2)/(q); 
                          var2 = eye(6) - 1/3*(I'*I) - (2/3).*(n*n');
                          var3 = (2*G)/(sqrt(3)).*n+K*obj.beta.*I';
                          var4 = (2*G)/(sqrt(3)).*n+obj.alpha*K.*I';
                          var5 = G+obj.alpha*obj.beta*K+obj.epsilon^2*obj.h;
                          D(:,:,i) = D(:,:,i) - var1*var2-var3*((var4)/(var5))';
                      end   
                  else
                      continue
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
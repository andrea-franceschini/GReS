classdef Elastic < handle
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
    function obj = Elastic(fID, matFileName)
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
      sigmaIn = lambdacorr(obj, sigmaIn, nptGauss);
      % MODIFICA SN
%       for i = 1 : nptGauss
%         sigmaOut(i,:) = sigmaIn(i,:) + epsilon(i,:)*D;
%         DAll(:,:,i) = D;
%       end
      sigmaOut = sigmaIn + epsilon*D;
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

    function sigmaIn = lambdacorr(obj, sigmaIn, nptGauss)
              for i = 1:nptGauss
                  q = sqrt(0.5*((sigmaIn(i, 1)-sigmaIn(i, 2))^2+(sigmaIn(i, 1)- ...
                      sigmaIn(i, 3))^2+(sigmaIn(i, 2)-sigmaIn(i, 3))^2)+ ...
                      sigmaIn(i, 4)^2+sigmaIn(i, 5)^2+sigmaIn(i, 6)^2); %calcolo di q
                  %k = 0; % inizialmente: dovr
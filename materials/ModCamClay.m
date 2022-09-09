classdef ModCamClay < handle
  % Modified Cam Clay material class

  properties (Access = private)
    sigma_c = [];            %Preconsolidation stress
    lambda = [];             %Isotropic compression index
    k = [];                  %Isotropic swelling index 
    nu_ur = [];              %Unloading-reloading Poisson ratio
    e0 = [];                 %Initial void ratio
    phi = [];                %Friction angle
    sigma_init = [];         %Initial stress vector
  end

  methods (Access = public)
    function obj = ModCamClay(inputString)
      obj.setMaterialParameters(inputString);
    end
    
    % Function calculating the stiffness matrix using the class properties
    function D = getStiffnessMatrix(obj, varargin)
        %Calculating mean effective stresse
        p_init = 1/3*(obj.sigma_init(1)+obj.sigma_init(2)+obj.sigma_init(3));
        a = (obj.sigma_init(1)-obj.sigma_init(2))^2;
        b = (obj.sigma_init(2)-obj.sigma_init(3))^2;
        c = (obj.sigma_init(3)-obj.sigma_init(1))^2;
        q_init = sqrt(1/2*(a+b+c+6*(obj.sigma_init(4)^2+obj.sigma_init(5)^2+obj.sigma_init(6)^2)));
        
        %Calculating principal effective stresses
        e = (obj.sigma_init(1)-p_init)*(obj.sigma_init(2)-p_init)*(obj.sigma_init(3)-p_init);
        f = -(obj.sigma_init(2)-p_init)*(obj.sigma_init(6))^2-(obj.sigma_init(1)-p_init)*(obj.sigma_init(5))^2-(obj.sigma_init(3)-p_init)*(obj.sigma_init(4))^2;
        g = 2*obj.sigma_init(4)*obj.sigma_init(5)*obj.sigma_init(6);
        J = e + f + g;
        teta = 1/3*asin(27/2*J/(q_init^3));
        sigma1 = p_init+2/3*q_init*sin(teta-2/3*pi);
        sigma2 = p_init+2/3*q_init*sin(teta);
        sigma3 = p_init+2/3*q_init*sin(teta+2/3*pi);

     
        % Calculating M (tangent of CSL) based on inizial mean effective
        % stress
        if sigma1 <= sigma2 && sigma2 == sigma3
            M = (6*sin(obj.phi))/(3-sin(obj.phi));
        elseif sigma1 == sigma2 && sigma2 <= sigma3
            M = (6*sin(obj.phi))/(3+sin(obj.phi));
        else
            M = sqrt(3)*sin(obj.phi);
        end
        
        %Calculating material stiffness matrix 
        D = zeros(6,6);
        K = ((1+obj.e0)*p_init)/obj.k;
        G = (9*(1-2*obj.nu_ur))/(2*(1+obj.nu_ur))*K;
        D(1,1) = K;
        D(2,2) = K;
        D(3,3) = K;
        D(1,2) = K;
        D(1,3) = K;
        D(2,1) = K;
        D(2,3) = K;
        D(3,1) = K;
        D(3,2) = K;
        D(4,4) = G;
        D(5,5) = G;
        D(6,6) = G;
        
        D = zeros(6,6);
        K = ((1+obj.e0)*p_init)/obj.k;
        G = (9*(1-2*obj.nu_ur))/(2*(1+obj.nu_ur))*K;
        D(1,1) = K;
        D(2,2) = K;
        D(3,3) = K;
        D(1,2) = K;
        D(1,3) = K;
        D(2,1) = K;
        D(2,3) = K;
        D(3,1) = K;
        D(3,2) = K;
        D(4,4) = G;
        D(5,5) = G;
        D(6,6) = G;
        
        
        
        
        
        
        
    end
  end

  
  
  methods (Access = private)
      % Function that set the material parameters coming from "data"
      % (Materials) inside the vector "params"
    function setMaterialParameters(obj, inputString)
      words = strsplit(inputString, ' ');
      params = zeros(length(words),1);
      k = 0;
      for i = 1 : length(words)
        if (length(words{i}) > 0)
          k = k + 1;
          params(k) = sscanf(words{i}, '%e');
        end
      end
      % Object properties are assigned with the same order used in the input
      % file
      obj.lambda = params(1);
      obj.k = params(2);
      obj.nu_ur = params(3);
      obj.e0 = params(4);
    end
  end

end

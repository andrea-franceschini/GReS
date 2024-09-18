classdef DKPlastic < handle
    % Drucker-Prager costitutive law implementation
    %Il materiale e elastico fino a quando F<0: quando F<0 allora vengono
    %calcolati G, K, n, e infine Δλ per poter correggere le tensioni. Viene
    %inserito un ciclo while per verificare che F sia effettivamente >0
    %dopo esser stato  corretto
    % Questo in un controllo if che verifica che q^tr sia maggiore di zero.
    % Infine vengono corrette le tensioni
    properties (Access = public)
    end

    methods (Access = public)
        function obj = DKPlastic(fID, matFileName)
            obj.readMaterialParameters(fID, matFileName);
        end
        
        function [status] = initializeStatus(obj, sigma)
          nptGauss = size(sigma,1);
          status = zeros(nptGauss,2);
        end

        function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status)
          nptGauss = size(sigmaIn,1); %numero di vertici dell'elemento
          D = zeros(6);
          D([1 8 15]) = 1-obj.nu;
          D([2 3 7 9 13 14]) = obj.nu;
          D([22 29 36]) = (1-2*obj.nu)/2;
          D = obj.E/((1+obj.nu)*(1-2*obj.nu))*D;
          sigmaIn = lambdacorr(obj, sigmaIn, nptGauss);
          sigmaOut = sigmaIn + epsilon*D;
          DAll = repmat(D,[1, 1, nptGauss]);
        end
        
        function m = getMFactor(obj)
          m = obj.M;
        end

        function cM = getRockCompressibility(obj)
          cM = obj.cM;
        end
        
        function sigmaIn = lambdacorr(obj, sigmaIn, nptGauss)
              for i = 1:nptGauss
                  q = sqrt(0.5*((sigmaIn(i, 1)-sigmaIn(i, 2))^2+(sigmaIn(i, 1)- ...
                      sigmaIn(i, 3))^2+(sigmaIn(i, 2)-sigmaIn(i, 3))^2)+ ...
                      sigmaIn(i, 4)^2+sigmaIn(i, 5)^2+sigmaIn(i, 6)^2); %calcolo di q
                  %k = 0; % inizialmente: dovrebbe essere obj.h*εp ma la
                  %deformazione plastica e zero
                  p = sum(sigmaIn(i, 1:3)); %calcolo di p
                  force = q/1.732 + obj.alpha*p-obj.epsilon*(obj.co+k);
                  while force > 0
                      sigma = sigmaIn;
                      G = obj.E/(2*(1+obj.nu));
                      K = obj.E/(3*(1+obj.nu));
                      %s = sigmaIn(i, 1:3)- p/3;
                      n = (1.5/q).*(sigma(i, 1:3)- p/3);
                      % Δλ
                      lambdac = max(0, force/(G+obj.alpha*obj.beta*K+(obj.epsilon^2)*obj.h)); 
                      p = p-lambdac*K*obj.beta;
                      if q-lambdac*1.732*G < 0
                          force = obj.alpha*p - obj.epsilon*(obj.co); %no shear force
                          sigma(i, 1:6) = (p-lambdac*K*obj.beta).*ones(1,6);
                      else
                          q = q-lambdac*1.732*G;
                          force = q/1.732 + obj.alpha*p-obj.epsilon*(obj.co+k); 
                          sigma(i, 1:6) = sigma(i, 1:6) - lambdac*((2*G/1.732).*n+K*obj.beta);
                      end
                  end
               end
        end
    end
    methods (Access = private)
        function readMaterialParameters(obj, fID, matFileName)
            tmpVec = readDataInLine(fID, matFileName, 2);
            obj.E = tmpVec(1);
            obj.nu = tmpVec(2);
            obj.psi = tmpVec(3);
            obj.phi = tmpVec(4);
            obj.co = tmpVec(5); %coesione
            obj.h = tmpVec(6); %hardening parameter
            obj.M = obj.nu/(1-obj.nu);
            obj.cM = (1+obj.nu)*(1-2*obj.nu)/(obj.E*(1-obj.nu));
            obj.alpha = (3*tan(obj.phi))/(sqrt(9+12*tan(obj.phi)^2)); %sti tre posso aggiungerli in getstifnessmatrix?
            obj.beta = (3*tan(obj.psi))/(sqrt(9+12*tan(obj.psi)^2));
            obj.epsilon = 3/(sqrt(9+12*tan(obj.phi)^2));
        end
    end    
end
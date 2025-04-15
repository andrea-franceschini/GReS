classdef HypoElastic < handle
    % HYPOELASTIC ISOTROPIC material class

    properties (Access = private)
        % Poisson's ratio
        nu
        % Coefficients a, b for virgin compressibility
        a
        b
        % Coefficients a1, b1 for unload/reload compressibility
        a1
        b1
        % Preconsolidation vertical stress
        szmin
        D1
        % M factor (sigmax/sigmaz)
        M
        %vector with updated preconsolidation vertical in all GP
        szvec = []
        % save simulation time step
        time = []

    end

    methods (Access = public)
        % Class constructor method
        function obj = HypoElastic(fID, matFileName)
            % Calling the function to set the object properties
            obj.readMaterialParameters(fID, matFileName);
            obj.computeConstPart();
        end

        function [status] = initializeStatus(obj, sigma)
            nptGauss = size(sigma,1);
            status = zeros(nptGauss,2);
        end

        % Material stiffness matrix calculation using the object properties
        %     function D = getStiffnessMatrix(obj, varargin)
        %       if (nargin ~= 2)
        %         error('Error in calling the HypoElastic/getStiffnessMatrix method - INPUT: (sz)');
        %       end
        %       % Stiffness matrix
        %       % varargin{1} is sz, i.e., the vertical stress
        %       cm = getCompressibility(obj, varargin{1});
        %       D = (1/cm)*obj.D1;
        %     end

        %NEW gestStiffnessMatrix

        function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, ~, status, el, t)
            %       if (nargin ~= 2)
            %         error('Error in calling the HypoElastic/getStiffnessMatrix method - INPUT: (sz)');
            %       end
            % Stiffness matrix
            % varargin{1} is sz, i.e., the vertical stress
            nptGauss = size(sigmaIn,1);
            DAll = zeros(6,6,nptGauss);
            sigmaOut = zeros(nptGauss,6);
            ez = epsilon(:,3);
            szIn = sigmaIn(:,3);
            dsAdim = epsilon*obj.D1;
            l1 = nptGauss*(el-1)+1 ;
            %inizialize obj.szvec (first D calculation) and obj.time
            if length(obj.szvec) < l1
                obj.szvec(l1:l1+nptGauss-1) = obj.szmin;
                obj.time(el) = -1; %initial value of time
            end
            szVec = obj.szvec(l1:l1+nptGauss-1);
            [szVec] = updatePrestress(obj, t, szIn, szVec, el);
            obj.szvec(l1:l1+nptGauss-1) = szVec;
            for i = 1:nptGauss
                if szIn(i)<= szVec(i)
                    az = obj.a;
                    bz = obj.b;
                else
                    az = obj.a1;
                    bz = obj.b1;
                end
                [D, sOut] = solverHypoElastic(obj, epsilon(i,:), szIn(i), ...
                    dsAdim(i,:), sigmaIn(i,:), az, bz);
                sigmaOut(i,:) = sOut;
                DAll(:,:,i) = D;
                %update szvec property for next step
                %if the new vertical stress is smaller than the
                %initial preconsolidation value, then the latter is updated
            end




            %cM = getRockCompressibility(obj, sigmaIn, epsilon, t); %medium value in the element
            %D = obj.D1.*reshape(1./cM,1,1,[]);
            %sigmaOut = sigmaIn + epsilon*D;
            %DAll = repmat(D,[1, 1, nptGauss]);
        end





        %OLD GETSTIFFNESSMATRIX
        %     function [DAll, sigmaOut, status] = getStiffnessMatrix(obj, sigmaIn, epsilon, dt, status)
        % %       if (nargin ~= 2)
        % %         error('Error in calling the HypoElastic/getStiffnessMatrix method - INPUT: (sz)');
        % %       end
        %       % Stiffness matrix
        %       % varargin{1} is sz, i.e., the vertical stress
        %       nptGauss = size(sigmaIn,1);
        %       cM = getRockCompressibility(obj, sigmaIn, epsilon, t); %medium value in the element
        %       D = obj.D1.*reshape(1./cM,1,1,[]);
        %       sigmaOut = sigmaIn + epsilon*D;
        %       DAll = repmat(D,[1, 1, nptGauss]);
        %     end

        % Method that returns the M factor
        function m = getMFactor(obj)
            m = obj.M;
        end

    end

    methods (Access = private)
        % Assigning material parameters (check also the Materials class)
        % to object properties
        function readMaterialParameters(obj, fID, matFileName)
            tmpVec = readDataInLine(fID, matFileName, 6);
            % Assign object properties
            obj.nu = tmpVec(1);
            obj.a = tmpVec(2);
            obj.b = tmpVec(3);
            obj.a1 = tmpVec(4);
            obj.b1 = tmpVec(5);
            obj.szmin = tmpVec(6);
            %
            % Compute the M factor
            obj.M = obj.nu/(1-obj.nu);
        end

        function computeConstPart(obj)
            % Stiffness matrix
            obj.D1 = zeros(6);
            obj.D1([1 8 15]) = 1;
            obj.D1([2 3 7 9 13 14]) = obj.nu/(1-obj.nu);
            obj.D1([22 29 36]) = (1-2*obj.nu)/(2*(1-obj.nu));
        end

        function [D, sigmaOut] = solverHypoElastic(obj, epsilon, szIn, dsAdim, sigmaIn, az, bz)
            %compute tangent stiffness matrix and updated stress state
            %in Gauss Point
            ez = epsilon(3);
            %first step, no vertical deformation
            %if  ez==0
            cm = az*abs(szIn)^bz;
            D = (1/cm)*obj.D1;
            sigmaOut = sigmaIn + epsilon*D;
            %         else
            %             szOut = (((bz+1)/az)*ez + abs(szIn)^(bz+1))^(1/(bz+1)); %absolute value
            %             phi = (szOut - abs(szIn))/dsAdim(3);
            %             sigmaOut = sigmaIn + dsAdim*phi;
            %             dsigma = sigmaOut - sigmaIn;
            %             cm = az*szOut^bz;
            %             bvec = -(1/cm)*obj.D1(3,:);
            %             cvec = -(1/dsAdim(3)^2)*obj.D1(3,:);
            %             B = dsAdim'*bvec;
            %             C = dsAdim'*cvec;
            %             D = (dsigma(3)/dsAdim(3))*obj.D1+...
            %                 (1/dsAdim(3))*B + dsigma(3)*C;
            %         end
        end



        function  [szVec] = updatePrestress(obj, t, szIn, szVec, el)
            if t == obj.time(el)
                return
            else
                log = szIn < szVec';
                szVec = log.*szIn + ~log.*szVec';
                obj.time(el) = t;
            end
        end




        % Compressibility calculation
        function cM = getRockCompressibility(obj, sigmaIn, epsilon)
            sz = mean(sigmaIn(:,3));
            if sz<obj.szmin
                % Loading path
                cM = (obj.a).*(abs(sz)).^(obj.b);
            else
                % Unloading/reloading path
                cM = (obj.a1).*(abs(sz)).^(obj.b1);
            end
        end
    end
end
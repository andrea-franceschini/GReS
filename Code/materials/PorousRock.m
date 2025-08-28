classdef PorousRock < handle
    % POROUS ROCK material class

    properties (Access = private)
        %General properties:
        KVec                 % Vector of permeabilities
        % (upper triangular part ordered column-wise)
        poro                 % Porosity
        biot                 % Biot coefficient
        %     alpha                % Rock compressibility (can be replaced by the oedometer test compressibility Cm)
        specGrav             % Specific gravity of rock
        % Swr                  % Residual saturation of water
        Sr=0.;             % Residual saturation
        Ss=1.;             % Maximum saturation
    end

    methods (Access = public)
        % Class constructor method
        function obj = PorousRock(fID, model, matFileName)
            % Calling the function to set the object properties
            obj.readMaterialParameters(model, fID, matFileName);
        end

        % Function to get material porosity
        function poro = getPorosity(obj)
            poro = obj.poro;
        end

        % Function to get material porosity
        % function Swr = getWaterResSat(obj)
        %     Swr = obj.Swr;
        % end

        function Ss = getMaxSaturation(obj)
            %GETSS Function to get the maximun saturation of the fluid.
            Ss = obj.Ss;
        end

        function Sr = getResidualSaturation(obj)
            %GETSS Function to get the residual saturation of the fluid.
            Sr = obj.Sr;
        end

        function specGrav = getSpecificGravity(obj)
            specGrav = obj.specGrav;
        end


        % Function to get material porosity
        function biotCoeff = getBiotCoefficient(obj)
            biotCoeff = obj.biot;
        end

        % Function to get material permeability as a 3x3 matrix
        function K = getPermMatrix(obj)
            if length(obj.KVec) == 1
                K = diag(obj.KVec*ones(3,1));
            elseif length(obj.KVec) == 3
                K = diag(obj.KVec);
            else
                K = [obj.KVec(1) obj.KVec(2) obj.KVec(3);
                    obj.KVec(2) obj.KVec(4) obj.KVec(5);
                    obj.KVec(3) obj.KVec(5) obj.KVec(6)];
            end
        end

        function K = getPermVector(obj) % Inspired by MRST
            if length(obj.KVec) == 1
                K = [obj.KVec, 0, 0, 0, obj.KVec, 0, 0, 0, obj.KVec];
            elseif length(obj.KVec) == 3
                K = [obj.KVec(1), 0, 0, 0, obj.KVec(2), 0, 0, 0, obj.KVec(3)];
            else
                K = [obj.KVec(1), obj.KVec(2), obj.KVec(3), obj.KVec(2), obj.KVec(4), ...
                    obj.KVec(5), obj.KVec(3), obj.KVec(5), obj.KVec(6)];
            end
        end

        % Function to get rock compressibility
        %     function a = getRockCompressibility(obj)
        %       a = obj.alpha;
        %     end
    end

    methods (Access = private)
        % Assigning material parameters (check also the Materials class)
        % to object properties
        function readMaterialParameters(obj, model, fID, matFileName)
            tmpVec =  sscanf(fgetl(fID), '%e');
            %tmpVec = readDataInLine(fID, matFileName, 3);
            obj.poro = tmpVec(1);
            obj.specGrav = tmpVec(2);
            if model.isFlow() && model.isPoromechanics()
                obj.biot = tmpVec(3);
            end
            KTmp = zeros(6,1);
            tmpVec = readDataInLine(fID, matFileName, 3);
            KTmp(1:3) = tmpVec;
            tmpVec = readDataInLine(fID, matFileName, 2);
            KTmp(4:5) = tmpVec;
            tmpVec = readDataInLine(fID, matFileName, 1);
            KTmp(6) = tmpVec;
            if model.isVariabSatFlow()
                tmpVec = readDataInLine(fID, matFileName, 2);
                obj.Sr = tmpVec(1);
                obj.Ss = tmpVec(2);
            end
            %       % Preliminary check on the number of rows in each material block
            %       % and the number of parameters
            %       nEntry = size(block,1);
            %       if nEntry ~= 5
            %         error('Wrong number of input rows in material %s',block(1));
            %       end
            %       KTmp = zeros(6,1);
            %       for i=2:5
            %         strParams = strsplit(block(i));
            %         nEntry = size(strParams,2);
            %         err = false;
            %         switch i
            %           case 3
            %             if nEntry ~= 3; err = true; end
            %           case {2 4}
            %             if nEntry ~= 2; err = true; end
            %           case 5
            %             if nEntry ~= 1; err = true; end
            %         end
            %         %
            %         if err
            %           error('Wrong number of input parameters in material %s, row %d',block(1),i);
            %         end
            %         %
            %         params = str2double(strParams);
            %         switch i
            %           case 2
            %             obj.poro = params(1);
            % %             obj.alpha = params(2);
            %             obj.specGrav = params(2);
            %           case 3
            %             KTmp([1 2 3]) = [params(1), params(2), params(3)];
            %           case 4
            %             KTmp([4 5]) = [params(1), params(2)];
            %           case 5
            %             KTmp(6) = params(1);
            %         end
            %       end
            if all(KTmp([2 3 5]) == 0)
                if all(KTmp([4 6]) == KTmp(1))
                    obj.KVec = KTmp(1);
                else
                    obj.KVec = [KTmp(1); KTmp(4); KTmp(6)];
                end
            else
                obj.KVec = KTmp;
            end
            clear KTmp
            %
            % K needs to be SPD. It is symmetric by construction but is it also
            % Positive Definite?
            K = getPermMatrix(obj);
            eigv = eig(K);
            if any(eigv < length(eigv)*eps(max(eigv)))
                % Tolerance chosen following the hint in:
                % https://it.mathworks.com/help/matlab/math/determine-whether-matrix-is-positive-definite.html#DetermineWhetherMatrixIsSPDExample-3
                error('The permeability matrix for material %s is not positive definite',matFileName);
            end
        end
    end
end
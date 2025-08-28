classdef VanGenuchten < handle
    %VANGENUCHTEN class to implement the richard's equation for the
    % retantion curve and the relative permeability.
    % 
    % Details.:
    % Water saturation as function of the pressure.
    % 
    % $$ S(p) = (1+(\beta p)^{n})^{-m} $$,
    %
    % where:
    % - $$ S(p) $$ is the effective water saturation as function of the pressure.
    % - $$ p $$ is the pressure head.
    % - $$ \beta,\ n,\ m$$ are empirical values to ajust the curve.
    % Implemented:
    %         Se = (1 + (b*p)^n).^(-m);
    %   d(Se)/dp = ((-m*n*p^(-1))*(b*p).^n * (1+(b*p)^n)^(-1-m));
    % d2(Se)/dp2 = ((m*n*p^(-2)).*(b*p).^n * (1+(b*p)^n)^(-2-m) ...
    %              * (n*(m*(b*p)^n-1) + (b*p)^n + 1));
    % 
    % Relative hydraulic conductivity by Mualem's Model:
    % $$ K_{r}(p) = (1+(\beta p)^{n})^{-\kappa m} ...
    %               [1- (\beta p)^{n m} (1+(\beta p)^{n})^{-m} ]^{2} $$
    % conditions:
    % -   0 < m < 1
    % -   m = 1-1/n
    % where:
    % - $$ K_{r}(p) $$ relative hydraulic conductivity.
    % - $$ p $$ is the pressure head.
    % - $$ \beta,\ n,\ m$$ are empirical values to ajust the curve.
    % - $$ \kappa $$ parameter to ajust the curve.
    % Implemented:
    %  Kr = (1+(b*p)^n)^(-k*m)* (1-((b*p)^(n*m))*((1+(b*p)^n)^(-m)))^2;
    % dKr = -(m*n/p)*(((b*p)^(m*n))-(((b*p)^n+1)^m)) ...
    %       * (((b*p)^n+1)^(-m*(k+2)-1)) *((k*(b*p)^n)*((b*p)^(m*n) ...
    %       - ((b*p)^n + 1)^m) - 2*((b*p)^(m*n)));
    % 
    % Relative hydraulic conductivity by Burdine's Model:
    % $$ K_{r}(p) = (1+(\beta p)^{n})^{-2m} ...
    %               [1- (\beta p)^{n m} (1+(\beta p)^{n})^{-m} ] $$
    % conditions:
    % -   0 < m < 1
    % -   m = 1-2/n
    % -   n > 2
    % where:
    % - $$ K_{r}(p) $$ relative hydraulic conductivity.
    % - $$ p $$ is the pressure head.
    % - $$ \beta,\ n,\ m$$ are empirical values to ajust the curve.
    % Implemented:
    %  Kr = (1 + (b*p)^n)^(-2.*m) * (1- ((b*p)^(n*m))*((1+(b*p)^n)^(-m)));
    % dKr = -(m*n*(p^(-1))) * ((1+(b*p)^n)^(-3.*m-1)) ...
    %       * (2*((b*p)^n)*((1+(b*p)^n)^m) + (b*p)^(m*n) - 2*(b*p)^(n*m+n));
    % 
    % Reference:
    % - An Introduction to Reservoir Simulation Using Matlab/Octave -
    % Knut Andrea, 2019
    % - A Closed-form Equation for Predicting the Hydraulic Conductivity of
    % Unsaturated Soils - M.Th.Van Genuchten, 1980 
    % - On Describing and predicting the hydraulic properties of
    % unsaturated soils - M.Th.Van Genuchten and D.R. Nielsen, 1985

    properties
        n;               % Empirical value - adimensional
        beta;            % Empirical value - need to be the same dimension as pressure.
        kappa;           % Empirical value para van Genuchten-Mualem permeability curve.
    end

    properties (Access = private)
        betaCor=false;   % Flag to indicated the necessity of correction for the beta.
        presCor=true;    % Flag to indicated the necessity of correction for the pressure.
        modelType;       % Flag to indicated the model type
        retantionCurve;  % Storage the retantion curve.
        relPermCurve;    % Storage the relative permability curve.
        nonNegPressure;  % Storage a flag for the non negative pressure.
    end

    methods (Access = public)
        function obj = VanGenuchten(fID,matFileName,varargin)
            %VanGenuchtenMualem Construct an instance of this class

            % If number of arguments is greater than 2, passing values.
            if (nargin>2) && (nargin<7)
                obj.passingData(varargin{:});
            else
                obj.readFromFile(fID,matFileName);
            end
        end

        function [Sw, dSw, ddSw] = computeSwAnddSw(obj,pres)
            % COMPUTESWANDDSW Method to compute the relative saturation and
            % it's derivatives
            
            % Compute the effective or normalized saturation.
            if strcmp(obj.modelType,'Tabular')
                [Sw, dSw, ddSw] = obj.retantionCurve.interpTable(pres);
            else
                p = obj.presCorrection(pres);
                [Sw, dSw, ddSw] = obj.computeSaturation(p);
                dSw = varCorrection(obj,dSw);
                ddSw = varCorrection(obj,ddSw);
            end
        end

        function [Kr, dKr] = computeRelativePermeability(obj,pres)
            %computeLwAnddLw Method to compute the relative permeability

            % Compute the relative permeability.
            switch obj.modelType
                case 'Mualem'
                    p = obj.presCorrection(pres);
                    [Kr, dKr] = obj.computeRelativePermeabilityMualem(p);
                    dKr = varCorrection(obj,dKr);
                case 'Burdine'
                    p = obj.presCorrection(pres);
                    [Kr, dKr] = obj.computeRelativePermeabilityBurdine(p);
                    dKr = varCorrection(obj,dKr);
                case 'Tabular'
                    [Kr, dKr, ~] = obj.relPermCurve.interpTable(pres);
                otherwise
            end
        end
    end

    methods
        function mm = getm(obj)
            %GETM - return the empirical value of m in the van genuchten
            % curves.
            if (obj.modelType)
                mm = 1-1/obj.n;
            else
                mm = 1-2/obj.n;
            end
        end

        function obj = betaCorrection (obj,specWeight)
            % BETACORRECTION Method to correct the beta parameter as function
            % of the fluid specific weight.
            if obj.betaCor
                obj.beta = obj.beta/specWeight;
                obj.betaCor = false;
            end
        end

        function flag = needBetaCor (obj)
            % NEEDBETACOR Return True or False to the necessity to correct
            % the beta parameter.
            flag = obj.betaCor;
        end
    end

    methods (Access = private)
        function readFromFile(obj,fID,matFileName)
            %readFromFile function to read from the file and constructed
            % the material class.

            modType = readToken(fID,matFileName);
            switch modType
                case 'Mualem'
                    tmpVec = readDataInLine(fID, matFileName, 3);

                    % Assign object properties
                    obj.n = tmpVec(1);
                    obj.beta = tmpVec(2);
                    obj.kappa = tmpVec(3);

                    obj.modelType = 'Mualem';
                case 'Burdine'
                    tmpVec = readDataInLine(fID, matFileName, 2);

                    % Assign object properties
                    % obj.modelType = false;
                    obj.n = tmpVec(1);
                    obj.beta = tmpVec(2);

                    obj.modelType = 'Burdine';
                case 'Tabular'
                    obj.modelType = 'Tabular';
                    obj.retantionCurve = TabularCurve(fID, matFileName);
                    obj.relPermCurve = TabularCurve(fID, matFileName);
                otherwise
                    % Assign object properties from a table
                    obj.readMaterialParametersFromTable(modType);
                    obj.betaCor = true;
                    obj.presCor = true;
                    obj.modelType = 'Mualem';
            end
        end

        function passingData(obj,varargin)
            %readFromFile function to read from the file and constructed
            % the material class.

            modType=varargin{1};
            switch modType
                case 'Mualem'
                    % Assign object properties
                    obj.modelType = 'Mualem';
                    obj.n = varargin{2};
                    obj.beta = varargin{3};
                    obj.kappa = varargin{4};
                case 'Burdine'
                    % Assign object properties
                    obj.modelType = 'Burdine';
                    obj.n = varargin{2};
                    obj.beta = varargin{3};
                case 'Tabular'
                    obj.modelType = 'Tabular';
                    obj.retantionCurve = TabularCurve(varargin{2}, varargin{3});
                    obj.relPermCurve = TabularCurve(varargin{2}, varargin{3});
                otherwise
                    % Assign object properties from a table
                    obj.modelType = 'Mualem';
                    obj.readMaterialParametersFromTable(modType);
                    obj.betaCor = true;
                    obj.presCor = true;                    
            end
        end

        function readMaterialParametersFromTable(obj,soil)
            %readMaterialParametersFromTable Assigning a table with some
            % types of soil and it's correspondent experimental values.
            names = [ "Clay"; "Clayey_loam"; "Loam"; "Loamy_sand"; ...
                "Sand"; "Sandy_clay"; "Sandy_clayey_loam"; ...
                "Sandy_loam"; "Silt"; "Silty_clay"; "Silty_clayey_loam";...
                "Silty_loam"];

            % Finding the position of the soil.
            pos = matches(names,soil);

            % Table of adimensional values for the experimental parameters.
            list_n = [0.098; 0.151; 0.168; 0.242; 0.502; 0.082; 0.124; ...
                0.161; 0.225; 0.121; 0.182; 0.221];
            alpha = [-1.825; -1.801; -1.954; -1.459; -1.453; -1.476; ...
                -1.676; -1.574; -2.182; -1.790; -2.076; -2.296];

            % Correcting the table.
            nn=list_n(pos);
            aa=alpha(pos);

            % Test if the soil is describe in the table.
            if not(isscalar(nn))
                error('Unknown soil type');
                return;
            end

            % Saving the information.
            obj.n = 10^(nn);
            obj.beta = 10^(aa+2);
            obj.kappa = 0.5;
        end

        function p = presCorrection (obj,pres)
            % PRESCORRECTION Method to correct the capilari pressure to be
            % positive and non-negative.
            if (obj.presCor)
                p = -pres;
                obj.nonNegPressure = p<=0;
                p(obj.nonNegPressure) = 1e-15;
                % p(p>1e4)=1e4;
            else
                p=pres;
            end
        end

        function [Sw, dSw, ddSw] = computeSaturation(obj,p)
            % COMPUTESATURATION compute the effective saturation as
            % function of the pressure.
            nn = obj.n;
            mm = obj.getm();
            b = obj.beta;

            % Solving the saturation and it's derivatives.
            Sw   = (1 + (b*p).^nn).^(-mm);
            dSw  = ((-mm*nn*p.^(-1)).*(b*p).^nn ...
                .* (1+(b*p).^nn).^(-1-mm));
            ddSw = ((mm*nn*p.^(-2)).*(b*p).^nn ...
                .* (1+(b*p).^nn).^(-2-mm) ...
                .*(nn*(mm*(b*p).^nn -1) ...
                + (b*p).^nn + 1));
        end

        function [Kr, dKr] = computeRelativePermeabilityMualem(obj,pres)
            %computeKrMualem compute the relative permeability by Mualem
            %Model.
            nn = obj.n;
            mm = obj.getm();
            b = obj.beta;
            k = obj.kappa;
            p = pres;
            
            % Solving the relative permeability and it's derivative.
            Kr = (1+(b*p).^nn).^(-k*mm).* (1- ...
                ((b*p).^(nn*mm)).*((1+(b*p).^nn).^(-mm))).^2;

            dKr = -(mm*nn./p).*(((b*p).^(mm*nn)) - (((b*p).^nn +1).^mm)) ...
                .* (((b*p).^nn +1).^(-mm*(k+2)-1)) ...
                .*( (k*(b*p).^nn).*((b*p).^(mm*nn) ...
                - ((b*p).^nn + 1).^mm) - 2*((b*p).^(mm*nn)));
        end

        function [Kr, dKr] = computeRelativePermeabilityBurdine(obj,pres)
            %computeKrBurdine compute the relative permeability by Burdine
            %Model.
            nn = obj.n;
            mm = obj.getm();
            b = obj.beta;
            p = pres;

            % Solving the relative permeability and it's derivative.
            Kr = (1 + (b*p).^nn).^(-2.*mm) .* ...
                (1- ((b*p).^(nn*mm)).*((1+(b*p).^nn).^(-mm)));
            dKr = -(mm*nn*(p.^(-1))) .* ((1+(b.*p).^nn).^(-3.*mm-1)) ...
                .* (2.*((b.*p).^nn).*((1+(b.*p).^nn).^mm) + ...
                (b*p).^(mm.*nn) - 2.*(b.*p).^(nn.*mm+nn));
        end

        function var = varCorrection(obj,var)
            % VARCORRECTION Method to correct the derivative of a variable.
            if (obj.presCor)
                var(obj.nonNegPressure) = 0.;
            end
        end
    end
end
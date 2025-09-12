classdef TabularCurve < handle

    properties (Access = private)
        tabW
        derivW
        derivW2
        nPoints
    end

    methods (Access = public)
        % Class constructor method
        function obj = TabularCurve(fID, matFileName)
            % Calling the function to set the object properties
            obj.readTabularCurve(fID, matFileName);
        end

        function varargout = interpTable(obj,x)
            % flip the sign of the pressure.
            x = -x;

            % Locate the intervals.
            xNeg = x <= 0;
            xOutL = x > 0 & x < obj.tabW(1,1);
            xOutR = x > obj.tabW(obj.nPoints,1);
            xIn = x >= obj.tabW(1,1) & x <= obj.tabW(obj.nPoints,1);

            % Linear extrapolation for the values outside the table, a greater the 0.
            if obj.tabW(1,1) == 0, dyL = 0;
            else, dyL = (obj.tabW(1,2) - 1)/obj.tabW(1,1);
            end
            dyR = (obj.tabW(obj.nPoints,2) - obj.tabW(obj.nPoints-1,2))/(obj.tabW(obj.nPoints,1) - obj.tabW(obj.nPoints-1,1));

            varargout{1} = zeros(length(x),1);
            varargout{1}(xNeg) = 1;
            varargout{1}(xOutL) = 1 + dyL*x(xOutL);
            varargout{1}(xOutR) = obj.tabW(obj.nPoints,2) + dyR*(obj.tabW(obj.nPoints,1) - x(xOutR));
            varargout{1}(xIn) = interp1(obj.tabW(:,1),obj.tabW(:,2),x(xIn),'linear');
            
            % Check for negative values
            idOutR = find(xOutR);
            yNeg = varargout{1}(xOutR) < 0;
            varargout{1}(idOutR(yNeg)) = 0;

            if nargout > 1 % compute first derivative
                % Locate the intervals.
                xNeg = x <= 0;
                xOutL = x > 0 & x < obj.derivW(1,1);
                xOutR = x > obj.derivW(obj.nPoints-1,1);
                xIn = x >= obj.derivW(1,1) & x <= obj.derivW(obj.nPoints-1,1);

                dyL = -obj.derivW(1,2)/obj.derivW(1,1);

                varargout{2}(xNeg) = 0;
                varargout{2}(xOutL) = dyL*x(xOutL);
                varargout{2}(xOutR) = dyR;
                varargout{2}(xIn) = -interp1(obj.derivW(:,1),obj.derivW(:,2),x(xIn),'linear');
            end
            if nargout > 2 % compute second derivative
                xNeg = x <= 0;
                xOutL = x > 0 & x < obj.derivW2(1,1);
                xOutR = x > obj.derivW2(obj.nPoints-2,1);
                xIn = x >= obj.derivW2(1,1) & x <= obj.derivW2(obj.nPoints-2,1);

                dyL = -obj.derivW2(1,2)/obj.derivW2(1,1);
                
                varargout{3} = zeros(length(x),1);
                varargout{3}(xNeg) = 0;
                varargout{3}(xOutL) = dyL*x(xOutL);
                varargout{3}(xOutR) = 0;
                varargout{3}(xIn) = interp1(obj.derivW2(:,1),obj.derivW2(:,2),x(xIn),'linear');
            end
        end
    end

    methods (Access = private)
        function readTabularCurve(obj, fID, matFileName)
            curveFname = readToken(fID, matFileName);
            obj.tabW = load(curveFname);
            obj.nPoints = size(obj.tabW,1);
            avgLen = max(diff(obj.tabW(:,1)));
            Len = obj.tabW(1,1);
            if (Len~=0 && Len > avgLen)
                fprintf(['The distance between the origin and the first point \n' ...
                    'of the curve is greater than the maximun distance between \n' ...
                    'the points, which can make convergence more difficult.\n']);
            end

            % first derivative
            obj.derivW = zeros(obj.nPoints-1,1);
            obj.derivW(:,1) = 0.5*(obj.tabW(1:obj.nPoints-1,1)+obj.tabW(2:obj.nPoints,1));
            obj.derivW(:,2) = diff(obj.tabW(:,2))./diff(obj.tabW(:,1));
            % second derivative
            obj.derivW2 = zeros(obj.nPoints-2,1);
            obj.derivW2(:,1) = obj.tabW(2:obj.nPoints-1,1);
            obj.derivW2(:,2) = diff(obj.derivW(:,2))./diff(obj.derivW(:,1));

            %h1 = diff(obj.tabW(1:obj.nPoints-1,1));
            %h2 = diff(obj.tabW(2:obj.nPoints,1));
            %obj.derivW2(:,2) = (h1.*obj.tabW(3:obj.nPoints,2) - (h1+h2).*obj.tabW(2:obj.nPoints-1,2) + h2.*obj.tabW(1:obj.nPoints-2,2))./ ...
            %                   (0.5*h1.*h2.*(h1+h2));
            %norm(diff(obj.derivW(:,2))./diff(obj.derivW(:,1)) - obj.derivW2(:,2))

            % because of -x
            obj.derivW(:,2) = -obj.derivW(:,2);
        end
    end
end
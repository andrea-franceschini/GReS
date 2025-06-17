classdef SolverStatistics < handle
    %SOLVERSTATISTICS Class to save all the information of relation of the
    % nonlinear solver.

    properties
        itMaxNR
        relTol
        absTol
        listTimeStep = []
        listIterByTimeStep = []
        listRelError = []
        listAbsError = []
        listPosBackStep = []
        listNumBackStep = []        
    end

    properties (Access = private)
        saveRelError    % flag to indicate to save the relative error.
        saveAbsError    % flag to indicate to save the absolute error.
        saveBackStep    % flag to indicate to save the back step information.
        itIsBackStep    % flag to indicate if the last iteration was a back step.
    end

    methods (Access = public)
        function obj = SolverStatistics(maxIt,tolRel,tolAbs,flags)
            %SOLVERSTATISTICS Construct an instance of this class
            obj.itMaxNR = maxIt;
            obj.relTol = tolRel;
            obj.absTol = tolAbs;
            obj.saveRelError = flags(1);
            obj.saveAbsError = flags(2);
            obj.saveBackStep = flags(3);
            obj.itIsBackStep = false;
        end

        function saveIt(obj,time,absError,relError)
            %SAVEIT Update the list of information of the nonlinear solver

            % Check necessity to save the information.
            if (obj.checkNActive())
                return
            end

            % Save the information about each time step used.
            obj.listTimeStep = [obj.listTimeStep;time];

            % flag to save the information of the number of the time step;
            flag = true;

            % Save the relative error information.
            if (obj.saveRelError)
                obj.listRelError = [obj.listRelError;relError];
                obj.listIterByTimeStep = [obj.listIterByTimeStep;length(relError)];
                flag = false;
            end

            % Save the absolute error information.
            if (obj.saveAbsError)
                obj.listAbsError = [obj.listAbsError;absError];
                if (flag)
                    obj.listIterByTimeStep = [obj.listIterByTimeStep;length(absError)];
                end
            end

            % Update the type of iteration.
            obj.itIsBackStep = false;
        end

        function saveBackIt(obj)
            %SAVEBACKIT Update the information of the number of back steps

            % Check necessity to save the information.
            if (obj.checkNActive() || ~obj.saveBackStep)
                return
            end

            % Saving the information.
            if (obj.itIsBackStep)
                nIt = length(obj.listPosBackStep);
                obj.listNumBackStep(nIt) = obj.listNumBackStep(nIt)+1; 
            else
                nIt = length(obj.listIterByTimeStep);
                obj.listPosBackStep = [obj.listPosBackStep;nIt];
                obj.listNumBackStep = [obj.listNumBackStep;1];
            end
            % Update the type of iteration.
            obj.itIsBackStep = true;
        end

        function [values, interval] = findIntervalValues(obj,it,type)
            %FINDINTERVALVALUES return the interval with the values saved.
            % type = 0 - Relative Error
            %        1 - Absolute Error

            % Check if the information is saved.
            values = [];
            interval = [];
            if (type)
                if ~obj.saveAbsError
                    return
                end
            else
                if ~obj.saveRelError
                    return
                end
            end

            % Check if the iteration is validy.
            niter = length(obj.listIterByTimeStep);
            if it>niter && it>0
                return
            end

            % Find the interval to plot.
            interval = zeros(2,1);
            if (it==1)
                interval(1) = 1;
            else
                interval(1) = sum(obj.listIterByTimeStep(1:it-1))+1;
            end
            interval(2) = sum(obj.listIterByTimeStep(1:it));

            if (type)
                values = obj.listAbsError(interval(1):interval(2));
            else
                values = obj.listRelError(interval(1):interval(2));
            end
        end

        function plotRelError(obj,it)
            %PLOTRELERROR plot the relative error convergency graphic for
            %the iteration given.

            % Check if the information is saved.
            if ~obj.saveRelError
                return
            end

            % Check if the iteration is validy.
            niter = length(obj.listIterByTimeStep);
            if it>niter && it>0
                return
            end

            % Find the interval to plot.
            [values, ~] = findIntervalValues(obj,it,false);

            % Plot the graphic.
            figure();
            semilogy(values,'black','LineWidth',2,'MarkerSize',14);
            hold on;
            semilogy([1;obj.listIterByTimeStep(it)],[obj.relTol;obj.relTol],'black--','LineWidth',2,'MarkerSize',14);
            xlabel('Number of Iteration')
            ylabel('Relative Error')
            xlim([1,obj.listIterByTimeStep(it)]);
            % ylim([obj.relTol,1e0]);
            % legend(tstr, 'Location', 'northeast')
            set(gca,'FontName','Liberation Serif','FontSize',16,'XGrid','on','YGrid','on')
            hold off;
        end

        function plotAbsError(obj,it)
            %PLOTABSERROR plot the absolute error convergency graphic for
            %the iteration given.

            % Check if the information is saved.
            if ~obj.saveAbsError
                return
            end

            % Check if the iteration is validy.
            niter = length(obj.listIterByTimeStep);
            if it>niter && it>0
                return
            end

            % Find the interval to plot.
            [values, ~] = findIntervalValues(obj,it,true);

            % Plot the graphic.
            figure();
            semilogy(values,'black','LineWidth',2,'MarkerSize',14);
            % semilogy(obj.listAbsError(posI:posF),'black','LineWidth',2,'MarkerSize',14);
            hold on;
            semilogy([1;obj.listIterByTimeStep(it)],[obj.absTol;obj.absTol],'black--','LineWidth',2,'MarkerSize',14);
            xlabel('Number of Iteration')
            ylabel('Absolute Error')
            xlim([1,obj.listIterByTimeStep(it)]);
            % ylim([obj.relTol,1e0]);
            % legend(tstr, 'Location', 'northeast')
            set(gca,'FontName','Liberation Serif','FontSize',16,'XGrid','on','YGrid','on')
            hold off;
        end

        function plotNIterTime(obj)
            %PLOTNITERTIME plot a graphic with the number of iteration
            % against the time simulated.

            % Check if the information is saved.
            if obj.checkNActive()
                return
            end

            % Plot the graphic.
            figure();
            plot(obj.listTimeStep,obj.listIterByTimeStep-1,'-','LineWidth',2,'MarkerSize',14);
            xlabel('Time')
            ylabel('Number of Iteration')
            % xlim([1,obj.listIterByTimeStep(it)]);
            % ylim([obj.relTol,1e0]);
            % legend(tstr, 'Location', 'northeast')
            set(gca,'FontName','Liberation Serif','FontSize',16,'XGrid','on','YGrid','on')
        end

        function plotNIter(obj)
            %PLOTNITER plot a graphic with the number of iteration against
            % the time step.

            % Check if the information is saved.
            if obj.checkNActive()
                return
            end

            % Plot the graphic.
            figure();
            plot(obj.listIterByTimeStep-1,'-','LineWidth',2,'MarkerSize',14);
            xlabel('Time Step')
            ylabel('Number of Iteration')
            % xlim([1,obj.listIterByTimeStep(it)]);
            % ylim([obj.relTol,1e0]);
            % legend(tstr, 'Location', 'northeast')
            set(gca,'FontName','Liberation Serif','FontSize',16,'XGrid','on','YGrid','on')
        end

        function plotBackSteps(obj)
            %PLOTBACKSTEPS plot the location of the back step in relation
            %of the time step

            % Check if the information is saved.
            if (obj.checkNActive() && obj.saveBackStep)
                return
            end

            listPos = linspace(0,length(obj.listTimeStep),length(obj.listTimeStep)+1);
            backsteps = zeros(1,length(obj.listTimeStep)+1);
            [~, pos] = ismember(obj.listPosBackStep, listPos);
            backsteps(pos) = obj.listNumBackStep;

            % Plot the graphic.
            figure();
            plot(listPos,backsteps,'-','LineWidth',2,'MarkerSize',14);
            xlabel('Time Step')
            ylabel('Number of Back Steps')
            set(gca,'FontName','Liberation Serif','FontSize',16,'XGrid','on','YGrid','on')
        end

        function plotBackStepsTime(obj)
            %PLOTBACKSTEPS plot the location of the back step in relation
            %of the time simulated.

            % Check if the information is saved.
            if (obj.checkNActive() && obj.saveBackStep)
                return
            end

            listPos = linspace(0,length(obj.listTimeStep),length(obj.listTimeStep)+1);
            backsteps = zeros(1,length(obj.listTimeStep)+1);
            [~, pos] = ismember(obj.listPosBackStep, listPos);
            backsteps(pos) = obj.listNumBackStep;

            % Plot the graphic.
            figure();
            plot([0; obj.listTimeStep],backsteps,'-','LineWidth',2,'MarkerSize',14);
            % plot(listPos,backsteps,'-','LineWidth',2,'MarkerSize',14);
            xlabel('Time Step')
            ylabel('Number of Back Steps')
            set(gca,'FontName','Liberation Serif','FontSize',16,'XGrid','on','YGrid','on')
        end
    end

    methods (Access = private)
        function flag = checkNActive(obj)
            %CHECKNACTIVE Check if the object is not being using to store
            %statistics information of the solver.
            if ((~obj.saveRelError) && (~obj.saveAbsError) )%&& (~obj.saveRelError))
                flag = true;
            else
                flag = false;
            end
        end
    end
end
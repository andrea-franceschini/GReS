classdef SolverStatistics < handle
    %SOLVERSTATISTICS Class to save all the information of relation of the
    % non linear solver.

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
            arguments
                maxIt {mustBePositive}
                tolRel double {mustBePositive}
                tolAbs double {mustBePositive}
                flags (3,1) logical
            end
            obj.itMaxNR = maxIt;
            obj.relTol = tolRel;
            obj.absTol = tolAbs;
            obj.saveRelError = flags(1);
            obj.saveAbsError = flags(2);
            obj.saveBackStep = flags(3);
            obj.itIsBackStep = false;
        end

        function saveIt(obj,time,relError,absError)
            %SAVEIT Update the list of information of the no-linear solver

            if (obj.checkActive())
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
            if (obj.checkActive() || obj.saveBackStep)
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


    end

    methods (Access = private)
        function flag = checkActive(obj)
            %CHECKACTIVE Check if the object is being using to store
            %statistics information of the solver.
            if ((~obj.saveRelError) && (~obj.saveAbsError) )%&& (~obj.saveRelError))
                flag = true;
            else
                flag = false;
            end
        end

    end
end
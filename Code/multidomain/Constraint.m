classdef (Abstract) Constraint < handle
  %CONSTRAINT Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

    interface
    coupledVariables
    

    
  end


  methods (Abstract)

    registerConstraint(obj,varargin)

    assembleConstraint(obj,varargin)

    updateState(obj,solution)

    % update the output structures for printing purposes
    [cellData,pointData] = writeVTK(obj,t);

    % write history to MAT-file
    writeSolution(obj,t,tID);

  end


  
  methods
    function obj = Constraint(inputArg1,inputArg2)
      %CONSTRAINT Construct an instance of this class
      %   Detailed explanation goes here
      obj.Property1 = inputArg1 + inputArg2;
    end
    
    function outputArg = method1(obj,inputArg)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      outputArg = obj.Property1 + inputArg;
    end
  end
end


classdef (Abstract) interfaceSolver < handle
  % Interface to create an interface between domains in GReS

  % Any interfaceSolver should implemented everything needed to couple a
  % variable field across two non conforming lower dimensional (2D) interfaces

  % The interface can belong to two different domains 

  % The fields property lists a set of single physical moduls that are
  % implemented

  % If the physicsSolver is coupled:
  % fields is a string array with all the single physics fields
  % J and rhs are cell variable with the same size as the number of fields
 

  properties (Abstract, GetAccess=public, SetAccess=public)
    J
    rhs
  end


  properties (GetAccess=public, SetAccess=protected)
    domain
  end

  methods
    function obj = interfaceSolver(domain, inputStruct)

      % domain:  handle to the Discretizer object storing all the
      % information of the model

      % read xml field that must have the same name of the class itself
      obj.readInput(inputStruct.(class(obj)));
    end
  end

  methods (Abstract)

    readInput(obj,input);
    setup(obj); 
    % register the existing state variables and DoF fields
    assembleSystem(obj); % compute the matrix and the rhs
    updateState(obj); % update the state variables 
    finalizeState(obj); % for post a processing variable
    [cellData,pointData] = buildPrintStruct(obj); % save outuputs to structures

  end
end



end
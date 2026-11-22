classdef Model < handle
  % GReS model manager
  
  properties

    domains
    interfaces

    J               % jacobian matrix
    res             % residual vector

    globalDoF = struct('dofs',[],'numbVars',[])       % store offset for global dof mapping 

  end
  
  methods

    function obj = Model(domains,interfaces)

      setModel(domains,interfaces)

    end


    function setModel(obj,domains,interfaces)

      nD = numel(domains);

      if nD == 0
        error("A model requires at least one Domain object")
      end
   
      obj.domains = domains;

      if nargin < 2
          obj.interfaces = [];
      else
        obj.interfaces = interfaces;
      end
    
      dofOffset = zeros(numel(obj.domains)+numel(obj.interfaces));
      varOffset = zeros(numel(obj.domains)+numel(obj.interfaces));

      for iD = 1:nD

        dofs = obj.domains(iD).dofm;
        dofOffset(iD) = dofs.getNumbDoF;
        varOffset(iD) = dofs.getNumberOfVariables(obj);
      end

      for iI = 1:numel(obj.interfaces)

        dofs = obj.interfaces(iI).dofm;
        dofOffset(nD+iI) = dofs.getNumbDoF;
        varOffset(nD+nI) = dofs.getNumberOfVariables(obj);

        if ~all(ismember(obj.interfaces.domains,obj.domains))
          error("Interface %s connects domains that are not in the model")
        end
      end

      if numel(obj.domains) > numel(unique(obj.domains))
        error('Multiple instance of the same Domain cannot be added to a Model')
      end
      if numel(obj.interfaces) > numel(unique(obj.interfaces))
        error('Multiple instance of the same Domain cannot be added to a Model')
      end

      obj.globalDoF.dofs = [1;cumsum(dofOffset(1:end-1))];
      obj.globalDoF.numbVars = [1;cumsum(varOffset(1:end-1))];



    end


    function assemble(obj)
    end

    function applyBoundaryConditions(obj)
    end

    function updateState(obj)
    end


    function printVTK(obj)
    end

    function setJacobian(obj,caller,J,rowVar,colVar)


      rowId = obj.getModelId(caller,rowVar);
      colId = obj.getModelId(caller,colVar);

      if ~all(size(obj.J{rowId,colId}) == size(J))
        error('Size of Jacobian block being added does not match existing block size')
      end

      obj.J{rowId,colId} = obj.J{rowId,colId} + J; 

    end

    function J = getJacobian(obj,caller,rowVar,colVar)

      if nargin == 1
        J = obj.J;
        return
      end

      if nargin == 2
        error('Incorrect number of input arguments')
      end

      if nargin == 3
        colVar = rowVar;
      end

      rowId = obj.getModelId(caller,rowVar);
      colId = obj.getModelId(caller,colVar);

      J = obj.J{rowId,colId};

    end


 function addJacobian(obj,caller,J,rowVar,colVar)

      if nargin < 4
        colVar = rowVar;
      end

      rowId = obj.getModelId(caller,rowVar);
      colId = obj.getModelId(caller,colVar);

      obj.J{rowId,colId} = obj.J{rowId,colId} + J;

    end





  
  end


  methods (Access = private)

    function idOut = getModelId(obj,input,varName)

      % return domain Id or global variable Id in the model

      switch class(input)
        case "Domain"
          idOut = find(obj.domains == input);

          if nargin > 2
            varId = obj.domains(idOut).dofm.getVariableId(varName);
            idOut = obj.globalDoF.numbVars(idOut) + varId;
          end
          
        case "Interface"
          id = find(obj.interfaces == input);
          idOut = numel(obj.domains) + id;

          if nargin > 2
            varId = obj.interfaces(id).dofm.getVariableId(varName);
            idOut = obj.globalDoF.numbVars(idOut) + varId;
          end
          
      end
    end
  end
end


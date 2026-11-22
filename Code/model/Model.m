classdef Model < handle
  % GReS model manager

  properties (SetAccess = private)
    domains
    interfaces
    J
    res
    dofOffsets
    fieldOffsets
    fieldSizes
  end

  methods
    function obj = Model(domains,interfaces)

      if nargin < 2
        interfaces = [];
      end
      obj.setModel(domains,interfaces);

    end


    function setModel(obj,domains,interfaces)

      if isempty(domains)
        error("A model requires at least one Domain object")
      end
      if nargin < 3
        interfaces = [];
      end

      obj.domains = domains;
      obj.interfaces = reshape(interfaces,1,[]);

      for i = 1:numel(domains)
        domains(i).model = obj;
      end

      for i = 1:numel(obj.interfaces)
        interf = obj.interfaces(i);
        if ~all(ismember(interf.domains,domains))
          error("An interface connects domains that are not in the model")
        end
        interf.model = obj;
      end

      obj.setDoFLayout();

    end


    function assembleSystem(obj,dt)

      obj.clearSystem();

      for i = 1:numel(obj.domains)
        obj.domains(i).assembleSystem(dt);
      end

      for i = 1:numel(obj.interfaces)
        obj.interfaces(i).assembleSystem();
      end

    end


    function updateState(obj,du)

      owners = obj.getOwners();
      offsets = obj.dofOffsets;

      for i = 1:numel(owners)
        owner = owners{i};
        increment = du(offsets(i)+1:offsets(i+1));
        owner.updateState(increment);
      end

    end


    function advanceState(obj)

      for owner = obj.getOwners()
        owner{1}.advanceState();
      end

    end


    function goBackState(obj)

      for owner = obj.getOwners()
        owner{1}.goBackState();
      end

    end


    function changed = updateConfiguration(obj)

      changed = false;
      for owner = obj.getOwners()
        ownerChanged = owner{1}.updateConfiguration();
        changed = changed || ownerChanged;
      end

    end


    function resetConfiguration(obj)

      for owner = obj.getOwners()
        owner{1}.resetConfiguration();
      end

    end


    function dofs = getGlobalDoFs(obj,owner,field)

      if nargin < 3
        id = obj.getOwnerId(owner);
        dofs = (obj.dofOffsets(id)+1:obj.dofOffsets(id+1))';
      else
        ids = obj.getFieldIds(owner,field);
        starts = [0,cumsum(obj.fieldSizes)];
        dofs = zeros(0,1);
        for id = ids
          dofs = [dofs; (starts(id)+1:starts(id+1))'];
        end
      end

    end


    function setJacobian(obj,owner,J,varargin)

      [rows,cols] = obj.getJacobianFields(owner,varargin{:});
      obj.J(rows,cols) = mat2cell(sparse(J),obj.fieldSizes(rows),obj.fieldSizes(cols));

    end


    function addJacobian(obj,owner,J,varargin)

      [rows,cols] = obj.getJacobianFields(owner,varargin{:});
      blocks = mat2cell(sparse(J),obj.fieldSizes(rows),obj.fieldSizes(cols));
      obj.J(rows,cols) = cellfun(@plus,obj.J(rows,cols),blocks,'UniformOutput',false);

    end


    function J = getJacobian(obj,owner,varargin)

      if nargin == 1
        J = obj.J;
      else
        [rows,cols] = obj.getJacobianFields(owner,varargin{:});
        J = obj.J(rows,cols);
      end

    end


    function setResidual(obj,owner,res,field)

      rows = obj.getFieldIds(owner,field);
      obj.res(rows) = mat2cell(sparse(res(:)),obj.fieldSizes(rows),1);

    end


    function res = getResidual(obj,owner,varargin)

      if nargin == 1
        res = obj.res;
      else
        res = cell2mat(obj.res(obj.getFieldIds(owner,varargin{:})));
      end

    end


    function applyDirichlet(obj,owner,field,dofs,isSymmetric)

      fieldId = obj.getFieldIds(owner,field);
      dofs = unique(dofs);

      for other = 1:numel(obj.fieldSizes)
        if ~isempty(obj.J{fieldId,other})
          obj.J{fieldId,other}(dofs,:) = 0;
        end
        if isSymmetric && ~isempty(obj.J{other,fieldId})
          obj.J{other,fieldId}(:,dofs) = 0;
        end
      end

      n = obj.fieldSizes(fieldId);
      obj.J{fieldId,fieldId} = obj.J{fieldId,fieldId} + ...
        sparse(dofs,dofs,ones(size(dofs)),n,n);
      obj.res{fieldId}(dofs) = 0;

    end


    function addResidualAtDoFs(obj,dofs,values)

      offsets = [0,cumsum(obj.fieldSizes)];
      for field = 1:numel(obj.fieldSizes)
        inField = dofs > offsets(field) & dofs <= offsets(field+1);
        if any(inField)
          local = dofs(inField) - offsets(field);
          obj.res{field} = obj.res{field} + accumarray(local(:),...
            reshape(values(inField),[],1),[obj.fieldSizes(field),1],@sum,0,true);
        end
      end

    end
  end


  methods (Access = private)

    function owners = getOwners(obj)
      % Get all domains and interfaces in a single cell array.

      owners = [num2cell(reshape(obj.domains,1,[])),num2cell(obj.interfaces)];

    end


    function setDoFLayout(obj)

      % prepare the Jacobian and residual cell array and the global dof
      % maps

      owners = obj.getOwners();
      nFields = cellfun(@(owner) owner.dofm.getNumberOfVariables(),owners);
      nDofs = cellfun(@(owner) owner.getNumbDoF(),owners);
      obj.dofOffsets = [0,cumsum(nDofs)];
      obj.fieldOffsets = [0,cumsum(nFields)];
      obj.fieldSizes = zeros(1,sum(nFields));

      for i = 1:numel(owners)
        for field = 1:nFields(i)
          obj.fieldSizes(obj.fieldOffsets(i)+field) = owners{i}.dofm.getNumbDoF(field);
        end
      end

      obj.clearSystem();

    end


    function clearSystem(obj)

      sizes = obj.fieldSizes;
      nFields = numel(sizes);
      obj.J = cell(nFields);
      obj.res = cell(nFields,1);

      for i = 1:nFields
        obj.res{i} = sparse(sizes(i),1);
        for j = 1:nFields
          obj.J{i,j} = sparse(sizes(i),sizes(j));
        end
      end

    end


    function id = getOwnerId(obj,owner)

      % get id of the owner in the model (domain or interface)

      if isa(owner,'Domain')
        id = find(obj.domains == owner,1);
      else
        id = numel(obj.domains) + find(obj.interfaces == owner,1);
      end

    end


    function ids = getFieldIds(obj,owner,fields)

      % get the location of the fields of the owner in the global system

      id = obj.getOwnerId(owner);
      if nargin < 3
        ids = obj.fieldOffsets(id)+1:obj.fieldOffsets(id+1);
      else
        if ~isnumeric(fields)
          fields = string(fields);
        end
        ids = arrayfun(@(field) owner.dofm.getVariableId(field),reshape(fields,1,[]));
        ids = ids + obj.fieldOffsets(id);
      end

    end


    function [rows,cols] = getJacobianFields(obj,owner,rowField,colField,colOwner)

      % return row and column of id of input fields in the system
      % separate owner for row and column fields are admitted


      if nargin < 3
        rows = obj.getFieldIds(owner);
      else
        rows = obj.getFieldIds(owner,rowField);
      end
      if nargin < 4
        cols = rows;
      elseif nargin < 5
        cols = obj.getFieldIds(owner,colField);
      else
        cols = obj.getFieldIds(colOwner,colField);
      end

    end
  end
end

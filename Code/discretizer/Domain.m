classdef Domain < handle

  % General model class
  % This class stores all the information related to the problem at hand
  % and store a map to the existing domainSolver and interfaceSolver



  properties (GetAccess=private, SetAccess=private)
    physicsSolvers        % db for physics solver in the model
  end

  % domain specific properties
  properties (GetAccess=public, SetAccess=private)
    dofm                % degrees of freedom manager
    bcs                 % boundary conditions
    outstate            % printing utilities
    materials           % materials
    % grid struct
    grid = struct('topology',[],'faces',[],'cells',[]) 
    elements            % db for finite elements
    faces               % db for finite volumes
    variables           % list of fields with the same order of the global system
    simparams
  end

  properties (GetAccess=public, SetAccess=public)
    state               % class holding any state variable in the domain
    stateOld            % class holding any state variable in the domain at the last converged time step

    % block jacobian and rhs defined in the domain
    J
    rhs
  end


  methods (Access = public)
    function obj = Domain(varargin)

      % the database for physicsSolvers in the domain
      obj.physicsSolvers = containers.Map('KeyType','char','ValueType','any');

      obj.setDomain(varargin{:})

    end


    function setDomain(obj,varargin)

      % create a domain from a key-value set of parameters or an xml file

      if isscalar(varargin)
        % input is a struct coming from an xml file
        obj.createDomainFromFile(varargin{1});
      else 
        % input is a list of key-value pair
        nIn = numel(varargin);
        assert(nIn == 8,"Error in setDomain(): " + ...
          "Incorrect key-value inputs. Keys must be:" + ...
          "'grid', 'materials', 'boundaryConditions', 'output'");
        fixedInput = struct(varargin{1:6});

        % read fixed input parameters
        obj.grid = fixedInput.grid;
        obj.materials = fixedInput.materials;
        obj.bcs = fixedInput.boundaryConditions;
        obj.outstate = fixedInput.output;

        % create the DoFManager
        obj.dofm = DoFManagerNew(obj.grid);

        % create the State object
        obj.state = State();

      end

    end

    function addPhysicsSolver(obj,solverName,varargin)

      % add a physic solver to an existing domain

      % ways to manually add a physics solver 
      % 1) a scalar struct coming from an xml file
      % 2) a key-value list of input, with additional solver specific input at the end

      if isscalar(varargin)
        % input is a xml file
        solverInput = varargin{1};
      else 
        % input is a list of key-value pair of solver inputs
        solverInput = struct(varargin{:});
      end

      % create solver
      pSolv = feval(solverName,obj);
      pSolv.registerSolver(solverInput);

      % add solver to the database
      obj.physicsSolvers(class(pSolv)) = pSolv;

      % redefine size of block jacobian and rhs for new added variables
      nVars = getNumberOfVariables(obj.dofm);
      obj.J = cell(nVars,nVars);
      obj.rhs = cell(nVars,1);

    end

    function out = getPhysicsSolver(obj,id)
      out = obj.physicsSolvers(id);
    end

    function out = getState(obj)
      out = obj.state;
    end

    function out = getStateOld(obj)
      out = obj.stateOld;
    end


    function J = getJacobian(obj,varargin)
      % GETJACOBIAN Return the system Jacobian matrix
      %
      % Usage:
      %   J = getJacobian(obj)
      %       Returns the full Jacobian matrix.
      %
      %   J = getJacobian(obj, fieldList)
      %       Returns only the Jacobian blocks corresponding to the
      %       specified fields. Assumes the same fields for both rows and columns.
      %
      %   J = getJacobian(obj, rowFields, colFields)
      %       Returns the Jacobian blocks corresponding to the specified
      %       fields for rows and columns separately.
      %
      % Inputs:
      %   fieldList  - string or cell array of field names (for both rows and columns)
      %   rowFields  - string or cell array of field names for rows
      %   colFields  - string or cell array of field names for columns
      %
      % Output:
      %   J          - the assembled Jacobian matrix
      %
      % Notes:
      %   - If no input fields are provided, the full Jacobian is returned.
      %   - Row and column fields must correspond to existing variables in the system.

      if nargin == 1
        J = obj.shrinkBlockMatrix(obj.J);
      elseif nargin == 2
        id = obj.dofm.getVariableId(varargin{1});
        J = obj.J{id,id};
      elseif nargin == 3
        idRow = obj.dofm.getVariableId(varargin{1});
        idCol = obj.dofm.getVariableId(varargin{2});
        J = obj.J{idRow,idCol};
      else
        error("Too many input arguments")
      end

    end

    function rhs = getRhs(obj,varargin)
      % GETRHS Return the right-hand side vector
      %
      % Usage:
      %   rhs = getRhs(obj)               - returns full RHS
      %   rhs = getRhs(obj, fieldList)    - returns only specified fields
      %
      % Inputs:
      %   fieldList - string or cell array of field names
      %
      % Output:
      %   rhs       - RHS vector (subset or full)
      %
      % Notes:
      %   Only one field list is allowed; multiple fields will be concatenated

      if nargin == 1
        rhs = obj.shrinkBlockMatrix(obj.rhs);
      elseif nargin == 2
        id = obj.dofm.getVariableId(varargin{1});
        rhs = obj.rhs{dofList};
      else
        error("Too many input arguments")
      end
    end


    function createDomainFromFile(obj,inputStruct)

      meshFile = getXMLData(inputStruct,[],"Geometry");
      mesh = Mesh();
      mesh.importMesh(meshFile);

      ng = getXMLData(inputStruct,1,"GaussPoints");

      obj.grid.topology = mesh;
      obj.grid.cells = Elements(mesh,ng);
      obj.grid.faces = Faces(mesh);

      obj.materials = Materials(inputStruct);
      obj.bcs = Boundaries(inputStruct,obj.grid);
      obj.outstate = OutState(mesh,inputStruct);

      % add physics solvers
      solvers = inputStruct.Solver;
      solvNames = fieldnames(solvers);

      % create the DoFManager
      obj.dofm = DoFManagerNew(obj.grid.topology);

      % create the State object
      obj.state = State();

      for s = 1:numel(solvNames)
        % define solver (the xml field must match the name of the solver
        % class)
        obj.addPhysicsSolver(solvNames{s},solvers.(solvNames{s}));
      end

    end







  end

  methods (Static)

    function C = shrinkBlockMatrix(C)
      % SHRINKMATRIX makes a sparse block cell array compatible with cell2mat.
      %
      % - Removes rows/columns that contain only empty cells.
      % - Replaces empty cells in remaining rows/columns with sparse matrices
      %   of the correct size for their block row and block column.

      [nR, nC] = size(C);

      % detect empty row and columns
      hasRow = false(nR,1);
      hasCol = false(1,nC);

      for i = 1:nR
        for j = 1:nC
          if ~isempty(C{i,j})
            hasRow(i) = true;
            hasCol(j) = true;
          end
        end
      end

      % remove empty rows/cols
      C = C(hasRow,hasCol);
      [nR,nC] = size(C);

      % determine block size
      rowSize = zeros(nR,1);
      for i = 1:nR
        for j = 1:nC
          if ~isempty(C{i,j})
            [rowSize(i),~] = size(C{i,j});
            break
          end
        end
      end

      % block column sizes (number of cols in each block column)
      colSize = zeros(1,nC);
      for j = 1:nC
        for i = 1:nR
          if ~isempty(C{i,j})
            [~,colSize(j)] = size(C{i,j});
            break
          end
        end
      end

      % replace empty blocks with empty sparse matrix
      for i = 1:nR
        for j = 1:nC
          if isempty(C{i,j})
            C{i,j} = sparse(rowSize(i), colSize(j));
          end
        end
      end
    end

  end
end

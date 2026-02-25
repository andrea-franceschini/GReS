classdef Materials < handle
  % MATERIAL - General material class

  properties (Access = public)
    solid     % container of all solid materials in the model
    fluid     % container of all fluid materials in the model
    matMap    % array map cell tag to corresponding material id
  end

  methods (Access = public)
    % Class constructor method
    function obj = Materials(input)

      obj.solid = cell([]);
      obj.fluid = cell([]);
      obj.matMap = [];
      
      if nargin == 0
        return
      end

      % Read any input material if provided
      obj.addMaterials(input);
    end

    function addMaterials(obj,input)

      input = readInput(input);

      % order matters: some solid PorousRock properties depend on the fluid
      if isfield(input,"Fluid")
        addFluid(obj,input.Fluid);
      end

      if isfield(input,"Solid")
        for i = 1:numel(input.Solid)
          addSolid(obj,input.Solid(i));
        end
      end

    end


    function addSolid(obj,varargin)

        matID = numel(obj.solid)+1;

        default = struct('cellTags',[],...
                         'name',string(strcat('mat_',num2str(matID))),...
                         'Constitutive',missing,...
                         'PorousRock',missing,...
                         'Curves',missing);

        input = readInput(default,varargin{:});

        % update the material map
        if length(obj.matMap) > max(input.cellTags)
          assert(~any(obj.matMap(cellTags)),...
            "Cannot assign material %s. Material have been already assigned to cell tag",input.name)
        end

        obj.matMap(input.cellTags) = matID;
        
        mat.name = input.name;

        % add solid material to the database
        obj.solid{matID} = mat;

        if ~ismissing(input.Constitutive)
          obj.addConstitutiveLaw(mat.name,input.Constitutive);
        end
        if ~ismissing(input.PorousRock)
          obj.addPorousRock(mat.name,input.PorousRock);
        end
        if ~ismissing(input.Curves)
          obj.addCapillaryCurves(mat.name,input.Curves);
        end
    end


    function addFluid(obj,varargin)

      obj.fluid{1} = Fluid(varargin{:}); 

    end


    function mat = getMaterial(obj,cellTag)
      % get material based on the cellTag using matMap

      mat = obj.solid{obj.matMap(cellTag)};

    end


    function materialNames = getMaterialNames(obj)

      db = [obj.solid{:}];
      materialNames = [db.name];

    end

    function matID = getMaterialIDFromName(obj,name)

      assert(isstring(name),"Material name must be a valid string");
      
      matNames = getMaterialNames(obj);
      matID = find(matNames==name);
      
      assert(isscalar(matID),"Multiple materials with name %s have been defined",matNames(matID(1)));

    end

    function mat = getMaterialFromName(obj,name)

      matID = getMaterialIDFromName(obj,name);
      mat = obj.solid{matID};

    end

    function fluidMat = getFluid(obj)
      % fluid material is always stored as the last one in the database
      fluidMat = obj.fluid;
    end

    function addConstitutiveLaw(obj,matID,varargin)
      % add a constitutive law to a material with specified name of id

      % add constitutive law
      if ~isnumeric(matID)
        matID = getMaterialIDFromName(obj,matID);
      end

      if nargin < 4
        % when called from addSolid()
        constLaw = fieldnames(varargin{1});
        if isstruct(varargin{1})
          input = varargin{1};
          input = input.(constLaw{:});
        end
        obj.solid{matID}.ConstLaw = feval(constLaw{1},input);
      else
        % when called externally from addSolid()
        constLaw = varargin{1};
        obj.solid{matID}.ConstLaw = feval(constLaw,varargin{2:end});
      end

    end

    function addPorousRock(obj,matName,varargin)

      if isempty(obj.fluid) 
        error("PorousRock properties can be added only if a fluid phase is present")
      end

      matID = obj.getMaterialIDFromName(matName);
      obj.solid{matID}.PorousRock = PorousRock(varargin{:});

    end

    function addCapillaryCurves(obj,matName,varargin)

      % this will be called from the PorousRock
      matID = getMaterialIDFromName(matName);
      obj.solid{matID}.PorousRock.Curves = VanGenuchten(varargin{:});

    end

    function [status] = initializeStatus(obj,cTag,sigma)
      mat = obj.getMaterial(cTag).ConstLaw;
      [status] = mat.initializeStatus(sigma);
    end

    function [D, sigma, status] = updateMaterial(obj, cTag, sigma, epsilon, dt, status, el, t)
      % constitutive update
      mat = obj.getMaterial(cTag).ConstLaw;
      [D, sigma, status] = mat.getStiffnessMatrix(sigma, epsilon, dt, status, el);
    end

    % Destructor
    % function delete(obj)
    %   remove(obj.db,keys(obj.db));
    % end
  end

  methods (Access = private)
    % Reading the input file by material blocks
    function readInputFile(obj, input)

      if ~isstruct(input)
        input = readstruct(input,AttributeSuffix="");
      end

      if isfield(input,"Materials")
        input = input.Materials;
      end

      if isfield(input,"fileName")
        assert(isscalar(fieldnames(input)),"FileName, " + ...
          " must be a unique parameter.");
        input = readstruct(input.fileName,AttributeSuffix="");
      end

      if isfield(input,"Fluid")
         fluid = Fluid(input.Fluid);
         %name = getXMLData(input.Fluid,[],"name");
         %fluid.setName(name);
      end

      nSolid = 0;
      if isfield(input,"Solid")
        nSolid = numel(input.Solid);
        cTags = [input.Solid.cellTags];
        if isstring(cTags)
          cTags = str2num(strjoin(cTags));
        end
        maxCellTag = max(cTags);
        obj.matMap = zeros(maxCellTag,1);
        for i = 1:nSolid
          %mat.name = getXMLData(input.Solid(i),[],"name");
          cellTags = getXMLData(input.Solid(i),[],"cellTags");
          if any(obj.matMap(cellTags))
            existingCellTags = cellTags(obj.matMap(cellTags)~=0);
            error("Multiple materials assigned to cellTags %s",...
              sprintf("%i ", existingCellTags));
          end
          obj.matMap(cellTags) = i;
          if isfield(input.Solid(i),"Constitutive")
            if ~ismissing(input.Solid(i).Constitutive)
              constLaws = fieldnames(input.Solid(i).Constitutive);
              % assumes that the XML field has the same name as the
              % constitutive law class
              mat.ConstLaw = feval(constLaws{1},input.Solid(i).Constitutive);
            end
          end
          if isfield(input.Solid(i),"PorousRock")
            if ~ismissing(input.Solid(i).PorousRock)
              mat.PorousRock = PorousRock(input.Solid(i).PorousRock);
            end
          end
          if isfield(input.Solid(i),"Curves")
            if ~ismissing(input.Solid(i).Curves)
              mat.Curves = VanGenuchten(input.Solid(i).Curves);
            end
          end
          obj.solid{i} = mat;
        end
      end

      if ~all(obj.matMap)
        % check if some cell tags do not have any material assigned
        missingCellTags = find(obj.matMap==0);
        error("Missing material for cellTags %s",...
          sprintf("%i ", missingCellTags));
      end

      if isfield(input,"Fluid")
        obj.db{nSolid+1} = fluid;
      end

    end


  end

end

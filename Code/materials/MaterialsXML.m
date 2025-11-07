classdef MaterialsXML < handle
  % MATERIAL - General material class

  properties (Access = public)
    % Creation of a Dictionary object (faster than map)
    %db = configureDictionary("double","struct");
    % configureDictionary is not supported before 2023b
    db
    matMap
  end

  methods (Access = public)
    % Class constructor method
    function obj = MaterialsXML(fileName)
      obj.db = containers.Map('KeyType','double','ValueType','any');
      % Calling the function to read input data from file
      obj.readInputFile(fileName);
    end

    % Get the material defined by matIdentifier and check if it is a
    % key of the Map db
    function mat = getMaterial(obj,cellID)
      %
      % The preliminary check whether matID is key of db has been commented
      % since it is highly expensive
      %       if (obj.db.isKey(matID))
      [matID,~] = find(obj.matMap == cellID);
      assert(isscalar(matID),['Zero or Multiple materials assigned to elements',...
        ' with cellTags %i'], cellID)
      mat = obj.db(matID);
      %       else
      %       Displaying error message if matIdentifier is not a key
      %       of the map
      %         error('Material %s not present', matID);
      %       end
    end

    function fluidMat = getFluid(obj)
      % fluid material is always stored as the last one in the database
      fluidMat = obj.db(max(keys(obj.db)));
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
    function delete(obj)
      remove(obj.db,keys(obj.db));
    end
  end

  methods (Access = private)
    % Reading the input file by material blocks
    function readInputFile(obj, fileName)

      input = readstruct(fileName,AttributeSuffix="");
      input = input.materials;

      nSolid = 0;
      if isfield(input,"Solid")
        nSolid = numel(input.solid);
        maxCellTag = max([input.solid.cellTags]);
        obj.matMap = zeros(maxCellTag,1);
        for i = 1:nSolid
          cellTags = getXMLData(input.solid(i),[],"cellTags");
          if any(obj.matMap(cellTags))
            existingCellTags = cellTags(obj.matMap(cellTags)~=0);
            error("Multiple materials assigned to cellTags %s",...
              sprintf("%i ", existingCellTags));
          end
          obj.matMap(cellTags) = i;
          if isfield(input,"Constitutive")
            constLaws = fieldnames(input.constitutive);
            % assumes that the XML field has the same name as the
            % constitutive law class
            mat.constLaw = feval(constLaws{1},input.constitutive);
          end
          if isfield(input,"PorousRock")
            mat.PorousRock = PorousRock(input.PorousRock);
          end
          if isfield(input,"Curves")
            mat.Curves = VanGenuchten(input.Curves);
          end
          obj.db(i) = mat;
        end
      end

      if ~all(obj.matMap)
        % check if some cell tags do not have any material assigned
        missingCellTags = find(obj.matMap==0);
        error("Missing material for cellTags %s",...
          sprintf("%i ", missingCellTags));
      end

      if isfield(input,"Fluid")
        obj.db(nSolid+1) = Fluid(input.fluid);
      end

      fclose(fID);
    end


  end

end
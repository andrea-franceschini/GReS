% % % % classdef MaterialN
% % % %     properties (Access = private)
% % % %         db1 containers.Map
% % % %         db2 containers.Map
% % % %     end
% % % %
% % % %     properties (Dependent)
% % % %         db3 containers.Map  % Derived: info + links to db1/db2 keys
% % % %     end
% % % %
% % % %     methods
% % % %         %% Constructor
% % % %         function obj = MaterialN()
% % % %             if nargin > 0
% % % %                 obj.db1 = containers.Map('KeyType','double','ValueType','any');
% % % %                 obj.db2 = containers.Map('KeyType','double','ValueType','any');
% % % %             end
% % % %         end
% % % %
% % % %         %% Getters
% % % %         function value = getDB1(obj)
% % % %             value = obj.db1;
% % % %         end
% % % %
% % % %         function value = getDB2(obj)
% % % %             value = obj.db2;
% % % %         end
% % % %
% % % %         %% Setters
% % % %         function obj = setDB1(obj, newDB1)
% % % %             obj.db1 = newDB1;
% % % %         end
% % % %
% % % %         function obj = setDB2(obj, newDB2)
% % % %             obj.db2 = newDB2;
% % % %         end
% % % %
% % % %         %% db3 Getter — builds map with links only
% % % %         function value = get.db3(obj)
% % % %             combinedMap = containers.Map;
% % % %
% % % %             % Define the keys in db3 (can be custom)
% % % %             % Example: keys common to both db1 and db2
% % % %             commonKeys = intersect(obj.db1.keys, obj.db2.keys);
% % % %
% % % %             for k = commonKeys
% % % %                 key = k{1};
% % % %
% % % %                 % Custom summary info
% % % %                 info = sprintf("Linked key '%s' found in both db1 and db2", key);
% % % %
% % % %                 combinedMap(key) = struct( ...
% % % %                     'info', info, ...
% % % %                     'keyInDB1', key, ...
% % % %                     'keyInDB2', key ...
% % % %                 );
% % % %             end
% % % %
% % % %             value = combinedMap;
% % % %         end
% % % %     end
% % % % end


classdef MaterialN
   % MaterialN - Store the material proprieties.
   %
   % Details.:
   % The class 
   % Class compose by 3 database and 1 array of keywords.
   % - one for solid, unique for each subdomain.
   % - n fluid's for the intared domain.
   % - one correlated database
   properties (Access = private)
      db containers.Map % material database
      matMapSolid       % Map for the key related to solid material
      matMapFluid       % Map for the key related to fluid material
      matMapCoMat       % Map for the key related to correleted materials
   end

   properties (Access = public)
      sdomainKey %array % store the keys database to the solid proprieties in the domain
   end

   methods
      % Constructor
      function obj = MaterialN(model,fListName)
         % Create the basic struct to store.
         obj.db = containers.Map('KeyType', 'char', 'ValueType', 'any');
         matMapSolid = zeros(100,100);
         matMapFluid = zeros(20,20);
         matMapCoMat = zeros(50,50);

         % Read the list of the material to read
         listFilesMat = readstruct(fListName,AttributeSuffix="");

         % Checks consistency and completeness for material data
         MaterialN.checkMatData(model,listFilesMat);

         % Store the material
         obj.readData(listFilesMat);
      end

      % Method to store data
      function obj = storeData(obj, key, value)
         if ischar(key) || isStringScalar(key)
            obj.db(char(key)) = value;
         else
            error('Key must be a character vector or string scalar.');
         end
      end

      % Optional: Method to retrieve data
      function value = getData(obj, key)
         if isKey(obj.db, key)
            value = obj.db(key);
         else
            error('Key "%s" not found in the database.', key);
         end
      end

      % Optional: Display all keys
      function displayKeys(obj)
         disp('Stored keys:');
         disp(keys(obj.db));
      end
   end

   methods (Access = private)
      function readData(obj,fdata)
         %READDATA - function to read the material information file in
         %xml and construct the class object.


         % Store the information for the Solid class.
         nsolmat = length(fdata.Solid);
         obj.solidMat = repmat(struct('matID', [], 'cellTag', []), nsolmat, 1);
         for mat=1:nsolmat
            if isnumeric(fdata.Solid(mat).cellTag)
               cellTags = fdata.Solid(mat).cellTag;
            else
               cellTags = sscanf(fdata.Solid(mat).cellTag, '%i,');
            end
            obj.solidMat(mat).matID = mat;
            obj.solidMat(mat).cellTag = cellTags;
            readXMLSolid(obj,fdata.Solid(mat).path,mat,cellTags);
         end

         % Make the map with with subdomain(cellTag) and it's material.
         obj.cellsolid = zeros(max(vertcat(obj.solidMat.cellTag)),1);
         for mat=1:nsolmat
            temp=obj.solidMat(mat).cellTag;
            for cells=1:length(temp)
               if obj.cellsolid(mat)==0
                  obj.cellsolid(mat)=temp(cells);
               else
                  if obj.cellsolid(mat)~=temp(cells)
                     error('More than one type of element at the same subdomain!');
                  end
               end
            end
         end

         % Check if fluid information is given and the necessity in the physics.
         obj.fluidMat = repmat(struct('matID', 1, 'cellTag', 1), 1, 1);
         fFlowNeed = model.isSinglePhaseFlow() || model.isVariabSatFlow();
         if fFlowNeed && isfield(fdata,'Fluid')
            if isnumeric(fdata.Solid(mat).cellTag)
               cellTags = fdata.Fluid(1).cellTag;
            else
               cellTags = sscanf(fdata.Fluid(1).cellTag, '%i,');
            end
            readXMLFluid(obj,fdata.Fluid(1).path,mat,cellTags)
         end

         % Make the map with with subdomain(cellTag) and it's material.
         obj.cellfluid = zeros(max(vertcat(obj.fluidMat.cellTag)),1);
         for mat=1:nsolmat
            temp=obj.fluidMat(mat).cellTag;
            for cells=1:length(temp)
               if obj.cellfluid(mat)==0
                  obj.cellfluid(mat)=temp(cells);
               else
                  if obj.cellfluid(mat)~=temp(cells)
                     error('More than one type of element at the same subdomain!');
                  end
               end
            end
         end
      end

   end



   methods (Static = true)
      

      % TODO -- Implement this class to validated the material data.
      function checkMatData(model,fdata)
         % CHECKMATDATA - Checks consistency and completeness for material
         % data.

         % Store the information for the Solid class.
         subregionsTags = regexp(fdata.Tags, '\w+', 'match');
      end

   end

end
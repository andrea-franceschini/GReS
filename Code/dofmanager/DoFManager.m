classdef DoFManager < handle
   % Degree of freedom manager
   % Map each entity to a corresponding degree of freedom in the solution
   % system
   % 2 ordering can be specified in input:
   % Field-based ordering (default)
   % Domain-based ordering 
   properties (Access = private)
      model
      cellTags
      domainFields
      tag2subDomain
      fieldList
   end
   properties (Access = public)
      numEntsField
      numEntsSubdomain
      fields
      ordering
      nComp
   end
   
   methods
      function obj = DoFManager(mesh,model,varargin)
         obj.model = model;
         obj.cellTags = mesh.cellTag;
         obj.fields = struct('field',[],'subdomainEnts',[],...
            'subID',[],'scheme',[],'entCount',[]);
         obj.tag2subDomain = zeros(mesh.nCellTag,1);
         % deal with variable imput
         switch nargin
            case 2 % no subdomains defined
               obj.ordering = 'field';
               obj.tag2subDomain(:) = 1; % only one subdomain
               availFields = obj.model.getAvailPhysics();
               for field = (string(availFields))'
                  % assign all available fields to one subdomain
                  addField(obj,mesh,1,field);
               end
               obj.numEntsField = zeros(1,numel(obj.fields));
               obj.numEntsSubdomain = 0;
            case 3 % no ordering specification
               obj.ordering = 'field';
               obj.readInputFile(mesh,varargin{1});
            case 4
               obj.ordering = varargin{1};
               obj.readInputFile(mesh,varargin{2});
            otherwise 
               error('Too many input arguments for class DoFManager');
         end
         % finilize DoF Manager construction
         obj.finalizeDoFManager(mesh);
      end
      
      function readInputFile(obj,mesh,fName)
         % Read DoF manager and update field informations
         fID = openReadOnlyFile(fName);
         l = 0;
         subID = 0;
         while ~feof(fID)
            line = readToken(fID,fName);
            l = l+1;
            if ~isempty(line) && ~strcmp(line(1),'%')
               assert(strcmp(line,'<SubDomain>'),['Error in DoFManager' ...
                  'input file %s at line %i of input file: unexpected ' ...
                  'instruction'])
               subID = subID+1;
               nextline = readToken(fID,fName);
               l = l+1;
               while ~strcmp(nextline,'</SubDomain>')
                  switch nextline
                     case '<CellTag>'
                        cTag = sscanf(fgetl(fID), '%i');
                        l = l+1;
                     case '<Field>'
                        fldList = convertCharsToStrings(split(strip((strtok(fgetl(fID),'%')))));
                        l = l+1;                
                     otherwise
                        if ~isempty(nextline)
                           error('Invalid block name in %s DoFmanager input file',fName);
                        else
                           error(['Error in %s input file: Invalid blank line ' ...
                              'in <Domain> block'],fName);
                        end
                  end
                  nextline = readToken(fID,fName);
                  l = l+1;
               end
               assert(all(obj.tag2subDomain(cTag)==0),['Error in subdomain' ...
                  ' %i: CellTag already defined'],subID);
               obj.tag2subDomain(cTag) = subID;
               % add field to DoFManager
               for field = fldList'
                  addField(obj,mesh,subID,field);
               end
            end
         end
         obj.numEntsField = zeros(1,numel(obj.fields));
         obj.numEntsSubdomain = zeros(1,subID);
         obj.nComp = zeros(1,numel(obj.fields));
      end


      function addField(obj,mesh,subIDs,field)
         % prepare fields structure from input file informations
         isField = false;
         if isscalar(obj.fields)
            if(isempty(obj.fields.field))
               nFields = 0;
               fl = false;
            else
               fl = true;
            end
         else
            fl = true;
         end
         if fl
            existingFields = string({obj.fields.field});
            nFields = numel(existingFields);
            isField = ismember(existingFields,field);
         end
         if any(isField)         
            % field aleady defined - add cell tag to list
            k = find(isField);
            assert(~ismember(obj.fields(k).subID,subIDs));
            obj.fields(k).subID = [obj.fields(k).subID subIDs]; 
         else % defining new field - initialize structure
            k = nFields+1;
            obj.fields(k).field = field;
            obj.fields(k).subID = subIDs;
            if isFEMBased(obj.model,translatePhysic(field))
               nEnts = mesh.nNodes;
               obj.fields(k).scheme = 'FEM';
            elseif isFVTPFABased(obj.model,translatePhysic(field))
               nEnts = mesh.nCells;
               obj.fields(k).scheme = 'FV';
            end
            obj.fields(k).isEntActive = false(nEnts,1);
         end
      end

      function finalizeDoFManager(obj,mesh)
         % updated DoF structure with entitiy list for each subdomain
         for i = 1:numel(obj.fields)
            nSub = numel(obj.fields(i).subID);
            obj.fields(i).entCount = zeros(1,nSub);
            obj.fields(i).subdomainEnts = cell(nSub,1);
            obj.nComp(i) = componentNumber(mesh,obj.fields(i).field);
            for j = 1:nSub
               subID = obj.fields(i).subID(j);
               cTags = find(obj.tag2subDomain == subID);
               switch obj.fields(i).scheme
                  case 'FEM'
                     ents = unique(mesh.cells(ismember(mesh.cellTag,cTags),:));
                  case 'FV'
                     ents = find(ismember(mesh.cellTag,cTags));
               end
               % store only entities not already assigned to a subdomain
               idActiveEnt = obj.fields(i).isEntActive(ents) == false;
               obj.fields(i).isEntActive(ents) = true; % activate entities
               obj.fields(i).subdomainEnts{j} = ents(idActiveEnt);
               obj.fields(i).entCount(j) = sum(idActiveEnt);
               obj.numEntsSubdomain(subID) = obj.numEntsSubdomain(subID) + obj.fields(i).entCount(j);
               obj.numEntsField(i) = obj.numEntsField(i) + obj.fields(i).entCount(j);
               obj.fieldList = string({obj.fields.field});
            end
         end
      end

      function dofs = getDoF(obj,field,varargin)
         % varargin: entity list, if empty all dofs are returned
         % return global DoF numbering based on entity ID and field
         assert(isscalar(string(field)),['Only one input field is allowed when ' ...
            'calling getDoF method']);
         fldId = getFieldId(obj,field);
         nc = obj.nComp(fldId);
         % get local numbering within the field
         switch obj.ordering
            case 'field'
               if isempty(varargin) 
                  % all dofs of input field
                  dofs = dofId((1:obj.numEntsField(fldId))',nc);
               else
                  dofs = getLocalDoF(obj,varargin{1},field);
               end
               nDoF = obj.nComp.*obj.numEntsField;
               dofs = dofs+sum(nDoF(1:fldId-1));
            case 'domain'
               error('Domain-based ordering of DoF not yet implemented')
         end
      end
      %
      function dofList = getLocalDoF(obj,entList,field)
         fldId = obj.getFieldId(field);
         nc = obj.nComp(fldId);
         % get local DoF numbering for entities within a field
         ents = getLocalEnts(obj,entList,field);
         dofList = dofId(ents,nc);
      end

      function dofList = getLocalEnts(obj,entList,field)
         % renumber entity id skipping inactive entities 
         fldId = obj.getFieldId(field);
         entList = reshape(entList,[],1);
         dofList = zeros(numel(entList),1);
         if ~all(obj.fields(fldId).isEntActive(entList))
            error('Inactive entity for field %s in input list',field);
         end
         % sorting entity list is way more efficient
         [entList,idSort] = sort(entList); 
         s = 0;
         k = 0;
         iPrev = 0;
         actEnt = obj.fields(fldId).isEntActive;
         for j = idSort'
            s = s+sum(actEnt(iPrev+1:entList(k+1)));
            iPrev = entList(k+1);
            dofList(j) = s;
            k = k+1;
         end
      end

      function fldDofs = getFieldDoF(obj,dofs,field)
         % recover active entitity indices from local dof numbering
         fldId = obj.getFieldId(field);
         actEnt = obj.fields(fldId).isEntActive;
         actDofs = find(actEnt);
         fldDofs = actDofs(dofs);
      end

      function activeSubs = getActiveSubdomain(obj,fieldList)
         % get subdomains where 1 or more fields are activated at the
         % same time
         % return [] if an input field is not available
         if ~all(ismember(fieldList,obj.fieldList))
             activeSubs = [];
             return
         end
         activeSubs = (1:max(obj.tag2subDomain))';
         for f = string(fieldList)
            fldId = getFieldId(obj,f);
            subField = obj.fields(fldId).subID;
            activeSubs = intersect(activeSubs,subField);
         end
      end

      function fldList = getFieldList(obj)
         fldList = obj.fieldList;
      end

      function fldId = getFieldId(obj,flds)
         % get field id associated to input fields
         fldId = find(strcmp(obj.fieldList,flds));
      end

      function cTags = getFieldCellTags(obj,field)
         % get cell ID where input field is active
         activeSubs = getActiveSubdomain(obj,field);
         cTags = find(ismember(obj.tag2subDomain,activeSubs));
      end
      
      function cellsID = getFieldCells(obj,field)
         cTags = getFieldCellTags(obj,field);
         cellsID = find(ismember(obj.cellTags,cTags));
      end

      function nComp = getDoFperEnt(obj,field)
         fldId = obj.getFieldId(field);
         nComp = obj.nComp(fldId);
      end

      function ents = getActiveEnts(obj,field)
         nc = obj.getDoFperEnt(field);
         id = obj.getFieldId(field);
         ents = dofId(find(obj.fields(id).isEntActive),nc);
      end

      function numDoF = getNumDoF(obj,field)
         fldId = obj.getFieldId(field);
         numEnts = obj.numEntsField(fldId);
         nc = obj.nComp(fldId);
         numDoF = nc*numEnts;
      end
   end

   methods (Access = private)
   end
end


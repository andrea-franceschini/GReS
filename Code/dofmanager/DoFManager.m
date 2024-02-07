classdef DoFManager < handle
    %DOFMANAGER Summary of this class goes here
    %   general subdomain manager
    %   associates a degree of freedom to nodes
    %test
    properties (Access = private)
        physicsList = "" %list of activated physics in the model
        ncomp %number of components for every physic in physicsList
        model 
        
    end

    properties (Access = public)
        subDomains %subDomain database 
        elemDofTable %table map elements ---> dof
        faceDofTable %table map faces ---> dof
        nodeDofTable %table nodes elements ---> dof
        numDof %array of total dofs for every subPhysic
        subPhysics %list of physics in the global solution vector
        subList %subDomain corresponding to every subPhysic
        subCells
    end
    
    methods
        function obj = DoFManager(mesh, model, varargin)
            obj.model = model;
            if isempty(varargin)
                fileName = [];
            else
                fileName = varargin{1};
            end
            obj.readInputFile(fileName,mesh); 
            %consider adding a default definition of subdomain if
            %fileName is empty (varargin)
            obj.buildDofTable(mesh)
            %obj.mergeSubDomains;
        end
        
        function readInputFile(obj,file,mesh)
            %building the subDomains structure and feature properties
            %in future some error messages will be added
            nSubMax = 10;
            obj.subDomains = repmat(struct('regions',[],'physics',[],'coupling',[]),nSubMax,1);
            if ~isempty(file)
                fID = obj.openReadOnlyFile(file);            
                assert(~feof(fID),'No SubDomain defined in the model');
                subID = 0;
                token  = [];
                while ~strcmp(token,'End')
                     subID = subID +1;
                     line = fgetl(fID);
                     domainRegions = sscanf(line, '%i');
                     physics = convertCharsToStrings(split(strip((strtok(fgetl(fID),'%')))));
                     if (~ismember(physics, ["Poromechanics","SPFlow","Thermal"]))
                         error(['Unknown physics in %s input file\n', ...
                         'Accepted physics for DoF Manager input file are: \n',...
                         'Poromechanics, SPFlow, Thermo'], file);
                     end
                     if physics(1) == "Coupled"
                         if length(physics) < 3
                             error('Coupled models need at least two physics')
                         else
                             obj.subDomains(subID).coupling = true;
                             physics = physics(2:end);
                         end
                     else
                         obj.subDomains(subID).coupling = false;
                     end
                     physics = rename(physics);
                     obj.subDomains(subID).regions = domainRegions;
                     obj.subDomains(subID).physics = physics;
                     tmp = ~ismember(physics,obj.physicsList);
                     obj.physicsList = [obj.physicsList; physics(tmp)];
                     token = fgetl(fID);
                end
                    obj.physicsList = obj.physicsList(2:end);
                    obj.subDomains =  obj.subDomains(1:subID); 
            else 
                    tmp = [isPoromechanics(obj.model); isFlow(obj.model)];
                    physics = ["Poro";"Flow"];
                    physics = physics(tmp);
                    obj.physicsList = physics;
                    obj.subDomains(1).physics = physics;
                    obj.subDomains(1).regions = (1:mesh.nCellTag)';                    
                    obj.subDomains = obj.subDomains(1);
            end
            obj.ncomp = ones(length(obj.physicsList),1)+(mesh.nDim-1)*ismember(obj.physicsList,"Poro");
            checkModelCompatibility(obj)
        end
            
        function buildDofTable(obj, mesh)
            %building matrices that map any entity (elements, faces,
            %nodes) with the corresponing DoF in the final system. 
            %code has to deal with "overlapping" subDomain physics
            %DofTable(:,:,1) ---> global dofs
            %DofTable(:,:,2) ---> local dofs
            %DofTbale(:,:,3) ----> physic index in global solution vector
            nPh = length(obj.physicsList)+ (mesh.nDim-1)*any(obj.physicsList == "Poro");
            %convert the following matrix in sparse format (after
            %testing)
            obj.subCells = sparse(mesh.nCells, length(obj.subDomains));
            obj.elemDofTable = zeros(mesh.nCells,nPh,3);
            obj.nodeDofTable = zeros(mesh.nNodes,nPh,3);
            numDofTable = zeros(length(obj.subDomains),length(obj.physicsList));
            dofCount = 0;
            phCount = 0;
            for subID = 1:length(obj.subDomains)
                obj.subCells(:,subID) = ismember(mesh.cellTag,obj.subDomains(subID).regions);
                for i=1:length(obj.physicsList)
                    %get column of dofTable taken by physics(i)
                    if any(obj.subDomains(subID).physics == obj.physicsList(i))
                        col = obj.getColTable(obj.physicsList(i));
                        if isFEMBased(obj.model,obj.physicsList(i)) %node centered discretization
                        %count number of zero entries for the nodesList for
                        %the specified physics in DofTable
                        nodesList = unique(mesh.cells(ismember(mesh.cellTag,obj.subDomains(subID).regions),:));
                        nodesList = nodesList(obj.nodeDofTable(nodesList, col(1))==0); %removing nodes with DOF already assigned
                        if ~isempty(nodesList)
                            numDofTable(subID,i) = length(nodesList)*obj.ncomp(i);
                            localDofs =  reshape(1:numDofTable(subID,i),obj.ncomp(i),length(nodesList));
                            globalDofs = dofCount + localDofs;
                            obj.nodeDofTable(nodesList,col,1) = globalDofs';
                            obj.nodeDofTable(nodesList,col,2) = localDofs';
                            phCount = phCount + 1;
                            obj.nodeDofTable(nodesList,col,3) = phCount;
                        end
                        elseif isFVTPFABased(obj.model,obj.physicsList(i))
                        elemList = find(ismember(mesh.cellTag,obj.subDomains(subID).regions));
                        if ~isempty(elemList)
                            numDofTable(subID,i) = length(elemList)*obj.ncomp(i);
                            localDofs = reshape(1:numDofTable(subID,i),obj.ncomp(i),length(elemList));
                            globalDofs = dofCount + localDofs; 
                            obj.elemDofTable(elemList,col,1) = globalDofs';
                            obj.elemDofTable(elemList,col,2) = localDofs';
                            phCount = phCount + 1;
                            obj.elemDofTable(elemList,col,3) = phCount;
                        end
                        end
                        dofCount = dofCount + numDofTable(subID,i);
                    end
                end
            end  
            [tmprow,tmpcol] = find(numDofTable' ~= 0);
            obj.subPhysics = obj.physicsList(tmprow);
            obj.subList = tmpcol;
            obj.numDof = nonzeros(numDofTable');
        end
        
        function dofs = getDoF(obj,physic,varargin)
            %Return Global DOFs associated to entities
            %Needed to map global Dofs with local solution arrays
            physic = translatePhysic(physic);
            col = obj.getColTable(physic);
            if isFEMBased(obj.model,physic) 
                if isempty(varargin)
                    dofs = (obj.nodeDofTable(:,col,1))';
                else
                    entities = varargin{1};
                    dofs = (obj.nodeDofTable(entities,col,1))';
                end
            elseif isFVTPFABased(obj.model,physic)
                if isempty(varargin)
                    dofs = (obj.elemDofTable(:,col,1))';
                else 
                    entities = varargin{1};
                    dofs = (obj.elemDofTable(entities,col,1))';
                end
            end
            dofs = dofs(:);
            dofs = nonzeros(dofs);
        end


        function subPh = getSubPhysic(obj,physic,varargin)
            %Return Physic ID associated to entities
            %Needed to map global Dofs with local solution arrays
            physic = translatePhysic(physic);
            col = obj.getColTable(physic);
            if isFEMBased(obj.model,physic) 
                if isempty(varargin)
                    subPh = (obj.nodeDofTable(:,col,3))';
                else
                    entities = varargin{1};
                    subPh = (obj.nodeDofTable(entities,col,3))';
                end
            elseif isFVTPFABased(obj.model,physic)
                if isempty(varargin)
                    subPh = (obj.elemDofTable(:,col,3))';
                else 
                    entities = varargin{1};
                    subPh = (obj.elemDofTable(entities,col,3))';
                end
            end
            subPh = subPh(:);
        end
        
        function dofs = getLocDoF(obj,physic,varargin)
            %Return Local subPhysics DOFs associated to entities
            %Needed for matrix assembly
            if ischar(physic)
                physic = translatePhysic(physic);
                col = obj.getColTable(physic);
                if isFEMBased(obj.model,physic) 
                    if isempty(varargin)
                        dofs = (obj.nodeDofTable(:,col,2))';
                    else
                        entities = varargin{1};
                        dofs = (obj.nodeDofTable(entities,col,2))';
                    end
                elseif isFVTPFABased(obj.model,physic)
                    if isempty(varargin)
                        dofs = (obj.elemDofTable(:,col,2))';
                    else 
                        entities = varargin{1};
                        dofs = (obj.elemDofTable(entities,col,2))';
                    end
                end
            else %if input is subPhysic ID
                if isFEMBased(obj.model,obj.subPhysics(physic))
                    [row,~] = find(obj.nodeDofTable(:,:,3)==physic);
                elseif isFVTPFABased(obj.model,obj.subPhysics(physic))
                    [row,~] = find(obj.elemDofTable(:,:,3)==physic);
                end
                    ncomps = obj.ncomp(obj.physicsList == obj.subPhysics(physic));
                    row = (reshape(row,[],ncomps))';
                    tmp = [0;1;2];
                    tmp = tmp(1:ncomps);
                    dofs = row(:)*length(tmp)-repmat(flip(tmp),size(row,2),1);     
            end
            dofs = dofs(:);
            dofs = nonzeros(dofs);
        end

        function ents = getEntities(obj,physic)
            %return indices of solution vector associated to input physic
            physic = translatePhysic(physic);
            col = obj.getColTable(physic);
            ncom = length(col);
            if isFEMBased(obj.model,physic)
                [ents,~] = find(obj.nodeDofTable(:,col) > 0);
            elseif isFVTPFABased(obj.model,physic)
                [ents,~] = find(obj.elemDofTable(:,col) > 0);
            end
            ents = (reshape(ents,[],ncom))';
            tmp = [0;1;2];
            tmp = tmp(1:ncom);
            ents = ents(:)*length(tmp)-repmat(flip(tmp),size(ents,2),1);     
        end

        function tags = getCellTag(obj,physic)
            tags = [];
            for i = 1:length(obj.subDomains)
                if any(strcmp(obj.subDomains(i).physics,physic))
                    tags = [tags; obj.subDomains.regions(:)];
                end
            end
        end

        function tags = getSubCells(obj,physic)
            tags = [];
            for i = 1:length(obj.subDomains)
                if any(strcmp(obj.subDomains(i).physics,physic))
                    tags = [tags; obj.subDomains.regions(:)];
                end
            end
        end
            
        
        function fID = openReadOnlyFile(obj,fName)
            fID = fopen(fName, 'r');
             if fID == -1
               error('File %s does not exist in the directory',fName);
             end
        end
        
        function col = getColTable(obj,physic)     
            i = find(obj.physicsList == physic);
            if i==1
                col = 1:obj.ncomp(i);   
            else
                col = (1:obj.ncomp(i))+sum(obj.ncomp(1:i-1));
            end
        end
                
         function table = getDofTables(obj)
             table.elemTable = obj.elemDofTable;
             table.nodeTable = obj.nodeDofTable;
         end
         
         function checkModelCompatibility(obj)
             if any(strcmp(obj.physicsList,"Poro"))
                 assert(isPoromechanics(obj.model),['Poromechanics physical ',...
                     'model not defined in ModelType class']);
              
             end
            
             if any(strcmp(obj.physicsList,"Flow"))
                 assert(isFlow(obj.model),['Flow physical ',...
                     'model not defined in ModelType class']);
             end    
         end


         function locList = glob2loc(obj, globList, ent)
             % map a global degree of freedome to its local counterpart
             switch ent
                 case 'node'
                    [row, col] = find(ismember(obj.nodeDofTable(:,:,1),globList'));
                    locTab = obj.nodeDofTable(:,:,2);
                    locList = locTab(sub2ind(size(locTab), row, col));
                 case 'elem'
                    [row, col] = find(ismember(obj.elemDofTable(:,:,1),globList'));
                    locTab = obj.elemDofTable(:,:,2);
                    locList = locTab(sub2ind(size(locTab), row, col));               
             end
         end

         function globList = loc2glob(obj, locList, ent)
             % map a local degree of freedome to its global counterpart
             switch ent
                 case 'node'
                    [row, col] = find(ismember(obj.nodeDofTable(:,:,2),locList'));
                    globTab = obj.nodeDofTable(:,:,1);
                    globList = globTab(sub2ind(size(globTab), row, col));
                 case 'elem'
                    [row, col] = find(ismember(obj.elemDofTable(:,:,2),locList'));
                    globTab = obj.elemDofTable(:,:,1);
                    globList = globTab(sub2ind(size(globTab), row, col));                  
             end
         end

         function subList = loc2sub(obj, locList, ent)
             % map a local dof to its physic ID
             switch ent
                 case 'node'
                    [row, col] = find(ismember(obj.nodeDofTable(:,:,2), locList'));
                    subTab = obj.nodeDofTable(:,:,3);
                    subList = subTab(sub2ind(size(subTab), row, col));
                 case 'elem'
                    [row, col] = find(ismember(obj.elemDofTable(:,:,2), locList'));
                    subTab = obj.elemDofTable(:,:,3);
                    subList = subTab(sub2ind(size(subTab), row, col));                
             end                         
         end

         function subList = glob2sub(obj, globList, ent)
             % map a local dof to its physic ID
             switch ent
                 case 'node'
                     [row, col] = find(ismember(obj.nodeDofTable(:,:,1), globList'));
                     subTab = obj.nodeDofTable(:,:,3);
                     subList = subTab(sub2ind(size(subTab), row, col));
                 case 'elem'
                     [row, col] = find(ismember(obj.elemDofTable(:,:,1), globList'));
                     subTab = obj.elemDofTable(:,:,3);
                     subList = subTab(sub2ind(size(subTab), row, col));
             end
         end
    end
end


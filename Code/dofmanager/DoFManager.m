classdef DoFManager < handle
    %DOFMANAGER Summary of this class goes here
    %   general subdomain manager
    %   associates a degree of freedom to nodes
    %test
    properties (Access = private)
        physicsList = ""    %list of activated physics in the model
        ncomp               %number of components for every physic in physicsList
        model 
        
    end

    properties (Access = public)
        subDomains      % subDomain database
        elem2dof        % table map elements ---> dof
        face2dof        % table map faces ---> dof
        nod2dof         % table nodes elements ---> dof
        numDof          % array of total dofs for every subPhysic
        subPhysics      % list of physics in the global solution vector
        subList         % subDomain corresponding to every subPhysic
        subCells        % list of cells of every subdomain
        glob2block     % map global dofs to field Tag
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
            obj.buildDofTable(mesh)
        end
        
        function readInputFile(obj,file,mesh)
            %building the subDomains structure and feature properties
            %in future some error messages will be added
            nSubMax = 10;
            obj.subDomains = repmat(struct('regions',[],'physics',[]),nSubMax,1);
            if ~isempty(file)
                fID = obj.openReadOnlyFile(file);            
                assert(~feof(fID),'No SubDomain defined in the model');
                subID = 0;
                token  = [];
                while ~strcmp(token,'End')
                     subID = subID +1;
                     line = fgetl(fID);
                     domainRegions = sscanf(line, '%i');
                     assert(all(ismember(domainRegions, mesh.cellTag)), ['Invalid' ...
                         'Region tag for subDomain %i'], subID)
                     physics = convertCharsToStrings(split(strip((strtok(fgetl(fID),'%')))));
                     if (~ismember(physics, ["Poromechanics","SPFlow","VSFlow"]))
                         error(['Unknown physics in %s input file\n', ...
                         'Accepted physics for DoF Manager input file are: \n',...
                         'Poromechanics, SPFlow, VSFlow'], file);
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
                    tmp = [isPoromechanics(obj.model); isSinglePhaseFlow(obj.model); isVariabSatFlow(obj.model)];
                    physics = ["Poro";"SPFlow";"VSFlow"];
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
            % row index --> entity tag
            % column index --> component tag
            % entry (~,1:end) --> global DoF
            nPh = length(obj.physicsList) + (mesh.nDim-1)*any(obj.physicsList == "Poro");
            %convert the following matrix in sparse format (after
            %testing)
            obj.subCells = sparse(mesh.nCells, length(obj.subDomains));
            obj.elem2dof = zeros(mesh.nCells,nPh);
            obj.nod2dof = zeros(mesh.nNodes,nPh);
            numDofTable = zeros(length(obj.subDomains),length(obj.physicsList));
            dofCount = 0;
            for subID = 1:length(obj.subDomains)
                obj.subCells(:,subID) = ismember(mesh.cellTag,obj.subDomains(subID).regions);
                for i=1:length(obj.physicsList)
                    %get column of dofTable taken by physics(i)
                    if any(obj.subDomains(subID).physics == obj.physicsList(i))
                        col = obj.getColTable(obj.physicsList(i));
                        if isFEMBased(obj.model,translatePhysic(obj.physicsList(i), obj.model)) %node centered discretization
                            %count number of zero entries for the nodesList for
                            %the specified physics in DofTable
                            nodesList = unique(mesh.cells(ismember(mesh.cellTag,obj.subDomains(subID).regions),:));
                            nodesList = nodesList(obj.nod2dof(nodesList, col(1))==0); %removing nodes with DOF already assigned
                            if ~isempty(nodesList)
                                numDofTable(subID,i) = length(nodesList)*obj.ncomp(i);
                                globalDofs = dofCount + reshape(1:numDofTable(subID,i),obj.ncomp(i),length(nodesList));
                                obj.nod2dof(nodesList,col) = globalDofs';
                            end
                        elseif isFVTPFABased(obj.model,translatePhysic(obj.physicsList(i),obj.model))
                            elemList = find(ismember(mesh.cellTag,obj.subDomains(subID).regions));
                            if ~isempty(elemList)
                                numDofTable(subID,i) = length(elemList)*obj.ncomp(i);
                                globalDofs = dofCount + reshape(1:numDofTable(subID,i),obj.ncomp(i),length(elemList));
                                obj.elem2dof(elemList,col) = globalDofs';
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
            obj.glob2block = repelem(1:length(obj.subPhysics), obj.numDof)';
        end
        
        function dofs = getDoF(obj,physic,varargin)
            %Return Global DOFs associated to entities
            %Needed to map global Dofs with local solution arrays
            col = obj.getColTable(physic);
            physic = translatePhysic(physic, obj.model);
            if isFEMBased(obj.model,physic) 
                if isempty(varargin)
                    dofs = (obj.nod2dof(:,col))';
                else
                    entities = varargin{1};
                    dofs = (obj.nod2dof(entities,col))';
                end
            elseif isFVTPFABased(obj.model,physic)
                if isempty(varargin)
                    dofs = (obj.elem2dof(:,col))';
                else 
                    entities = varargin{1};
                    dofs = (obj.elem2dof(entities,col))';
                end
            end
            dofs = dofs(:);
            dofs = nonzeros(dofs);
        end

        function fldDofs = ent2field(obj, fld, ents)
            % INPUT: entity ID, physic
            % OUTPUT: field dof numbering
            glob = obj.getDoF(fld, ents);
            k = cumsum(obj.numDof.*(~strcmp(obj.subPhysics, fld)));
            blkID = obj.glob2block(glob);
            tmp = [0; k(1:end-1)];
            fldDofs = glob - tmp(blkID);
        end

        function dofs = field2block(obj,fld)
            % INPUT: field ID
            % OUTPUT: block dof numbering corresponding to input field
            tmp = strcmp(obj.subPhysics, obj.subPhysics(fld));
            s = zeros(length(obj.numDof),1);
            s(tmp) = cumsum(obj.numDof(tmp));
            s = [0; s];
            dofs = ((s(fld)+1):(s(fld)+obj.numDof(fld)))';
        end


        function dofs = glob2blockDoF(obj, list)
            list = list(:); % make sure is a column vector
            % Return block numbering from global dof numbering
            s = cumsum(obj.numDof);
            s = [0; s(1:end-1)];
            dofs = list - s(obj.glob2block(list));
        end

        % function dofs = field2dof(obj,field)
        %     dofs = find(ismember(obj.glob2block,field));
        % end

        function ents = field2ent(obj, field)
            % INPUT: field string
            % OUTPUT: state array indices
            col = obj.getColTable(field);
            field = translatePhysic(field, obj.model);
            ncom = length(col);
            if isFEMBased(obj.model, field)
                ents = find(obj.nod2dof(:,col(1))); % non ordered active ents
                [~, i] = sort(nonzeros(obj.nod2dof(:,col(1))));
                ents = ents(i);
            elseif isFVTPFABased(obj.model, field)
                ents = find(obj.elem2dof(:,col(1))); % non ordered ents
                [~, i] = sort(nonzeros(obj.elem2dof(:,col(1))));
                ents = ents(i);
            end
            tmp = [0;1;2];
            tmp = tmp(1:ncom);
            ents = repelem(ents(:)*ncom,ncom)-repmat(flip(tmp),length(unique(ents)),1);
        end
        
        function ents = getEntities(obj, physic, varargin)
            %return indices of solution vector associated to input physic
            % or field (varargin)
            col = obj.getColTable(translatePhysic(physic, obj.model));
            ncom = length(col);
            if isFEMBased(obj.model,physic)
                [ents,~] = find(obj.nod2dof(:,col) > 0);
            elseif isFVTPFABased(obj.model,physic)
                [ents,~] = find(obj.elem2dof(:,col) > 0);
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

        function numbDof = getNumbDof(obj,physic)
            % get total number of active dofs for each physic
            id = strcmp(obj.subPhysics,physic);
            numbDof = sum(obj.numDof(id));
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
             table.elemTable = obj.elem2dof;
             table.nodeTable = obj.nod2dof;
         end
         
         function checkModelCompatibility(obj)
             for i = 1:length(obj.physicsList)
                 switch obj.physicsList(i)
                     case "Poro"
                         assert(isPoromechanics(obj.model),['Poromechanics physical ',...
                             'model not defined in ModelType class']);

                     case "SPFlow"
                         assert(isSinglePhaseFlow(obj.model),['Single Phase Flow physical ',...
                             'model not defined in ModelType class']);

                     case "VSFlow"
                         assert(isVariabSatFlow(obj.model),['Variably Saturated Flow physical ',...
                             'model not defined in ModelType class']);
                 end
             end
         end

    end
end


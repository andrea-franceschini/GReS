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
                    if length(physics) == 2
                        obj.subDomains(1).coupling = true;  
                    else
                        obj.subDomains(1).coupling = false;
                    end
                    
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
            obj.elemDofTable = zeros(mesh.nCells,nPh,3);
            obj.nodeDofTable = zeros(mesh.nNodes,nPh,3);
            numDofTable = zeros(length(obj.subDomains),length(obj.physicsList));
            dofCount = 0;
            phCount = 0;
            for subID = 1:length(obj.subDomains)
                for i=1:length(obj.physicsList)
                    phCount = phCount + 1;
                    %get column of dofTable taken by physics(i)
                    if any(obj.subDomains(subID).physics == obj.physicsList(i))
                        col = obj.getColTable(obj.physicsList(i));
                        if isFEMBased(obj.model,obj.physicsList(i)) %node centered discretization
                        %count number of zero entries for the nodesList for
                        %the specified physics in DofTable
                        nodesList = unique(mesh.cells(ismember(mesh.cellTag,obj.subDomains(subID).regions),:));
                        nodesList = nodesList(obj.nodeDofTable(nodesList, col(1))==0); %removing nodes with DOF already assigned
                        numDofTable(subID,i) = length(nodesList)*obj.ncomp(i);
                        localDofs =  reshape(1:numDofTable(subID,i),obj.ncomp(i),length(nodesList));
                        globalDofs = dofCount + localDofs;
                        obj.nodeDofTable(nodesList,col,1) = globalDofs';
                        obj.nodeDofTable(nodesList,col,2) = localDofs';
                        obj.nodeDofTable(nodesList,col,3) = phCount;
                        elseif isFVTPFABased(obj.model,obj.physicsList(i))
                        elemList = find(ismember(mesh.cellTag,obj.subDomains(subID).regions));
                        numDofTable(subID,i) = length(elemList)*obj.ncomp(i);
                        localDofs = reshape(1:numDofTable(subID,i),obj.ncomp(i),length(elemList));
                        globalDofs = dofCount + localDofs; 
                        obj.elemDofTable(elemList,col,1) = globalDofs';
                        obj.elemDofTable(elemList,col,2) = localDofs';
                        obj.elemDofTable(elemList,col,3) = phCount;
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
        end
        
        function dofs = getLocDoF(obj,physic,varargin)
            %Return Local subPhysics DOFs associated to entities
            %Needed for ApplyBCAndForces and matrix assembly
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
            dofs = dofs(:);
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
        
    
    
            
    end
end


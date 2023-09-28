classdef DoFManager < handle
    %DOFMANAGER Summary of this class goes here
    %   general subdomain manager
    %   associates a degree of freedom to nodes
    
    properties (Access = private)
        physicsList = ""
        ncomp
        model
    end
    
        properties (Access = public)
        subDomains
        elemDofTable
        faceDofTable
        nodeDofTable
        numDofTable
    end
    
    methods
        function obj = DoFManager(mesh,fileName,model)
            obj.model = model;
            obj.readInputFile(fileName,mesh);
            obj.buildDofTable(mesh)
            %obj.mergeSubDomains;
        end
        
        function readInputFile(obj,file,mesh)
            %building the subDomains structure and feature properties
            %in future some error messages will be added
            fID = obj.openReadOnlyFile(file);
            nSubMax = 10;
            obj.subDomains = repmat(struct('regions',[],'physics',[]),nSubMax,1);
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
             obj.subDomains(subID).regions = domainRegions;
             obj.subDomains(subID).physics = physics;
             tmp = ~ismember(physics,obj.physicsList);
             obj.physicsList = [obj.physicsList; physics(tmp)];
             token = fgetl(fID);
            end
            obj.physicsList = obj.physicsList(2:end);
            obj.ncomp = ones(length(obj.physicsList),1)+(mesh.nDim-1)*ismember(obj.physicsList,"Poromechanics");
            obj.subDomains =  obj.subDomains(1:subID);     
        end
            
        function buildDofTable(obj, mesh)
            %building matrices that map any entity (elements, faces,
            %nodes) with the corresponing DoF in the final system. 
            %code has to deal with "overlapping" subDomain physics
            nPh = length(obj.physicsList)+ (mesh.nDim-1)*any(obj.physicsList == "Poromechanics");
            %convert the following matrix in sparse format (after
            %testing)
            obj.elemDofTable = zeros(mesh.nCells,nPh);
            obj.nodeDofTable = zeros(mesh.nNodes,nPh);
            obj.numDofTable = zeros(length(obj.subDomains),length(obj.physicsList));
            dofCount = 0;
            for subID = 1:length(obj.subDomains)
                for i=1:length(obj.physicsList)
                    %get column of dofTable taken by physics(i)
                    if any(obj.subDomains(subID).physics == obj.physicsList(i))
                        if i==1
                            col = 1:obj.ncomp(i);   
                        else
                            col = (1:obj.ncomp(i))+sum(obj.ncomp(1:i-1));
                        end
                        if isFEMBased(obj.model,translatePhysic(obj.physicsList(i))) %node centered discretization
                        %count number of zero entries for the nodesList for
                        %the specified physics in DofTable
                        nodesList = unique(mesh.cells(ismember(mesh.cellTag,obj.subDomains(subID).regions),:));
                        nodesList = nodesList(obj.nodeDofTable(nodesList, col(1))==0); %removing nodes with DOF already assigned
                        obj.numDofTable(subID,i) = length(nodesList)*obj.ncomp(i);
                        dofs = dofCount + reshape(1:obj.numDofTable(subID,i),obj.ncomp(i),length(nodesList));
                        obj.nodeDofTable(nodesList, col) = dofs';
                        elseif isFVTPFABased(obj.model,translatePhysic(obj.physicsList(i)))
                        elemList = find(mesh.cellTag(ismember(mesh.cellTag,obj.subDomains(subID).regions),:));
                        obj.numDofTable(subID,i) = length(elemList)*obj.ncomp(i);
                        dofs = dofCount + reshape(1:obj.numDofTable(subID,i),obj.ncomp(i),length(elemList)); 
                        obj.elemDofTable(elemList, col) = dofs';
                        end
                        dofCount = dofCount + obj.numDofTable(subID,i);
                    end
                end
            end               
        end
        
        function dofs = getDoF(obj,entities,physic)
            phID = find(obj.physicsList == translatePhysic(physic));
            if phID == 1
                col = 1:obj.ncomp(phID);   
                else
                col = (1:obj.ncomp(phID))+sum(obj.ncomp(1:phID-1));
            end
            if isFEMBased(obj.model,physic)
                dofs = (obj.nodeDofTable(entities,col))';
            elseif isFVTPFABased(obj.model,physic)
                dofs = (obj.elemDofTable(entities,col))';
            end
            dofs = dofs(:);
        end
            
        
        function fID = openReadOnlyFile(obj,fName)
            fID = fopen(fName, 'r');
             if fID == -1
               error('File %s does not exist in the directory',fName);
             end
        end
                
         function table = getDofTables(obj)
             table.elemTable = obj.elemDofTable;
             table.nodeTable = obj.nodeDofTable;
         end
        
       
        
        
        end
    
            
end


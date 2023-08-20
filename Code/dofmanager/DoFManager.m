classdef DoFManager < handle
    %DOFMANAGER Summary of this class goes here
    %   general subdomain manager
    %   associates a degree of freedom to nodes
    
    properties
        db = containers.Map('KeyType','double','ValueType','any');
        nDoFs
    end
    
    methods
        function obj = DoFManager(mesh,fileName)
            obj.readInputFile(mesh,fileName);
            obj.mergeSubDomains;
        end
        
        function readInputFile(obj,mesh,file)
            fID = obj.openReadOnlyFile(file);
            assert(~feof(fID),'No SubDomain defined in the model');
            subID = 0;
            token  = [];
            while ~strcmp(token,'End')
             subID = subID +1;
             line = fgetl(fID);
             domainRegions = sscanf(line, '%i');
             physics = convertCharsToStrings(split(strip((strtok(fgetl(fID),'%')))));
             obj.db(subID) = SubDomain(mesh,domainRegions,physics);       
             token = fgetl(fID);
            end
           end
            
        
        
        
        function fID = openReadOnlyFile(obj,fName)
            fID = fopen(fName, 'r');
             if fID == -1
               error('File %s does not exist in the directory',fName);
             end
        end
        
        function sub = getSubDomain(obj,subID)
            sub = obj.db(subID);
        end
        
        function nodes = getSubNodes(obj,subID)
            nodes = obj.getSubDomain(subID).nodes;
        end
        
         function physics = getSubPhysics(obj,subID)
            physics = obj.getSubDomain(subID).physics;
         end
        
         function table = getSubTable(obj,subID)
             table = obj.getSubDomain(subID).dofTable;
         end
        
        
        
        function mergeSubDomains(obj)
            %manage DOFS shared by multiple subdomains   
            nSub = length(obj.db);
            for i=1:nSub-1
                ph1 = getSubPhysics(obj,i);
                n1 = getSubNodes(obj,i);
                for j=2:nSub
                    ph2 = getSubPhysics(obj,j);
                    n2 = getSubNodes(obj,j);
                    shared = intersect(n1,n2);
                    obj.getSubDomain(j).dofTable(ismember(n2,shared),ismember(ph1,ph2)) = 0; 
                end
            end
            
                    
                
            
        end
        
        
        end
    
            
    end


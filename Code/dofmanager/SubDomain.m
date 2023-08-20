classdef SubDomain < handle
    %SubDomain class. 
    %A subdomain is a container of dofs. 
    %if subdomains are defined, local matrices loops are performed 
    %within subdomains elements, not the entire domain.
    
    properties
        cells
        nodes
        physics
        mapDof
        dofTable
    end
    
    methods
        function obj = SubDomain(mesh,tags,physics)
            %SUBDOMAIN Construct an instance of this class
            %   Detailed explanation goes here
            obj.physics = physics;
            obj.cells = find(ismember(mesh.cellTag,tags));
            obj.nodes = unique(mesh.cells(obj.cells,:));
            obj.mapDof = physicsToDof(obj);
            obj.dofTable = computeDofTable(obj);
            
            
        end
        
        function mapDof = physicsToDof(obj)
            mapDof = zeros(length(obj.physics),1);
            for i=1:length(obj.mapDof)
                if obj.physics(i)=="Poromechanics"
                    mapDof(i) = 3;
                else
                    mapDof(i) = 1;
                end
            end
        end
    
        function table = computeDofTable(obj)
            table = ones(length(obj.nodes),length(obj.physics));
            table(:,obj.physics=="Poromechanics") = 3; 
        end
    end
end


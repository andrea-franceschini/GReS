classdef DofMap   
    properties
        nodeMapMaster
        nodeMapSlave
        nNMaster
        nNSlave
        nNIntMaster
        nNIntSlave
    end
    
    methods
        function obj = DofMap(meshMaster,nodesMaster,meshSlave,nodesSlave)
            % build dofMap
            % each node of a domain is mapped to its position in the
            % solution vector
            nodeMap = zeros(meshMaster.nNodes,1);
            obj.nNIntMaster = numel(nodesMaster);
            obj.nNMaster = numel(nodeMap)- obj.nNIntMaster;
            nodeMap(~ismember(1:meshMaster.nNodes,nodesMaster)) = (1:obj.nNMaster)';
            nodeMap(nodesMaster) =  (1:obj.nNIntMaster)';
            obj.nodeMapMaster = nodeMap;
            nodeMap = zeros(meshSlave.nNodes,1);
            obj.nNIntSlave = numel(nodesSlave);
            obj.nNSlave = numel(nodeMap)- obj.nNIntSlave;
            nodeMap(~ismember(1:meshSlave.nNodes,nodesSlave)) =  (1:obj.nNSlave)'+obj.nNMaster;
            nodeMap(nodesSlave) =  (1:obj.nNIntSlave)';
            obj.nodeMapSlave = nodeMap;
            obj.nodeMapMaster(nodesMaster) = obj.nodeMapMaster(nodesMaster)+...
                obj.nNMaster+obj.nNSlave;
            obj.nodeMapSlave(nodesSlave) = obj.nodeMapSlave(nodesSlave)+...
                obj.nNMaster+obj.nNSlave+obj.nNIntMaster;
        end

        function dofList = getDoF(obj,list,side,varargin)
            switch side
                case 'master'
                    nodList = obj.nodeMapMaster(list);
                case 'slave'
                    nodList = obj.nodeMapSlave(list);
                case 'lagrange'
                    nodList = obj.nodeMapSlave(list)+obj.nNIntSlave;
            end
            if isempty(varargin)
                dofList = ([2*nodList-1 2*nodList])';
                dofList = dofList(:);
            elseif strcmp(varargin{1},'x')
                dofList = 2*nodList-1;
            elseif strcmp(varargin{1},'y')
                dofList = 2*nodList;
            end
        end
    end

    methods (Static)
       function dofList = getCompDoF(list)
          % map dof numbering to component numbering
          dofList = zeros(2*numel(list),1);
          dofList(1:2:end) = 2*list-1;
          dofList(2:2:end) = 2*list;
       end
    end
end


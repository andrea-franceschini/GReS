classdef RBF_interpolation < handle
    % Perform Radial basis functions interpolation between different grids
    
    properties
        meshIn
        meshOut
        fiMM
        fiNM
        maxRad
        radiusList
        circumRadius
    end
    
    methods
        function obj = RBF_interpolation(mshIn, mshOut, nLinks)
            obj.meshIn = mshIn;      % Origin discretization
            obj.meshOut = mshOut;    % Destination discretization
            obj.computeFiMM(nLinks); % Interpolation matrix (return local radius if nLinks = 1) 
            obj.computeCrossMat();   
            obj.computeCircumRadius();   % maximum edge length of each cell 
        end

        function fOut = interpolate(obj,fIn)
            weightF = obj.fiMM\fIn;
            weightUnit = obj.fiMM\ones(obj.meshIn.nNodes,1);
            fOut = (obj.fiNM*weightF)./(obj.fiNM*weightUnit);
        end


        function P = compute_cross_mat(obj)
            % compute mortar method cross matrix on Slave interface only
            % make use of radial basis functions
        end

        function [wf, w1] = computeLocalInterpolant(obj, nId)
            % compute the weights for the RL-RBF interpolant around a
            % specific node of the grid
            rad = obj.radiusList(nId); % radius of local support
            
        end

        function obj = computeCircumRadius(obj)
            if any(obj.meshIn.coordinates(:,3)~=0)
                h = zeros(obj.meshIn.nCells,1);
                for el = 1:obj.meshIn.nCells
                    nodes = obj.meshIn.cells(el,:);
                    h(el) = max(vecnorm(obj.meshIn.coordinates(nodes,:) - obj.meshIn.coordinates([nodes(2:end) nodes(1)],:),2,2));
                end
            else
                h = zeros(obj.meshIn.nSurfaces,1);
                for el = 1:obj.meshIn.nSurfaces
                    nodes = obj.meshIn.surfaces(el,:);
                    h(el) = max(vecnorm(obj.meshIn.coordinates(nodes,:) - obj.meshIn.coordinates([nodes(2:end) nodes(1)],:),2,2));
                end
            end
            obj.circumRadius = h/2;
        end
        

        function  obj = computeFiMM(obj, radVal)
            % return fiMM matrix based on Wendland C2 RBF
            % radVal: number of links of neighboring nodes
            % compute distance table between nodes
            % loop trough nodes
            iiVec = []; jjVec = []; rbfVec = [];
            radList = zeros(obj.meshIn.nNodes,1);
            for i = 1:obj.meshIn.nNodes
                % compute local radius
                list = i;
                for count = 1:radVal
                    if any(obj.meshIn.coordinates(:,3)~=0)
                        cellNeigh = find(any(ismember(obj.meshIn.cells, list),2));
                        % get list of new neighboring nodes
                        ltmp = obj.meshIn.cells(cellNeigh,:);
                        list = unique([list; ltmp(:)]);
                    else
                        % find surfaces containing list of nodes
                        surfNeigh = find(any(ismember(obj.meshIn.surfaces, list),2));
                        % get list of new neighboring nodes
                        ltmp = obj.meshIn.surfaces(surfNeigh,:);
                        list = unique([list; ltmp(:)]);
                    end
                end
                tmp = obj.meshIn.coordinates(list,:) - obj.meshIn.coordinates(i,:);
                dist = sqrt(sum(tmp.^2,2));
                rad = max(dist);
                rbf = (1-dist./rad).^4.*(1+4*dist./rad);
                iiVec = [iiVec; repmat(i,length(list),1)];
                jjVec = [jjVec; list];
                rbfVec = [rbfVec; rbf];
                radList(i) = rad;
            end
            obj.fiMM = sparse(iiVec,jjVec,rbfVec,obj.meshIn.nNodes,obj.meshIn.nNodes);
            obj.maxRad = max(radList);
            obj.radiusList = radList;
        end

        function  obj = computeCrossMat(obj)
            % return fiNM matrix based on Wendland C2 RBF
            iiVec = []; jjVec = []; rbfVec = [];

            for i = 1:obj.meshOut.nNodes
                % for each node compute distances with all the others
                dist = obj.meshIn.coordinates - obj.meshOut.coordinates(i,:);
                dist = sqrt(sum(dist.^2,2));
                list = find(dist < obj.maxRad);
                dist = dist(list);
                rbf = (1-dist./obj.maxRad).^4.*(1+4*dist./obj.maxRad);
                iiVec = [iiVec; repmat(i,length(list),1)];
                jjVec = [jjVec; list];
                rbfVec = [rbfVec; rbf];
            end
            obj.fiNM = sparse(iiVec,jjVec,rbfVec,obj.meshOut.nNodes,obj.meshIn.nNodes);
        end


    end
end


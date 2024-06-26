classdef RayTracing< handle
    % Compute intersection between a ray and a surface
    % ray is described by origin and direction
    
    properties
        surf
        surfCentroids
        BBtree
        treeNodes
        leaf2elem
    end
    
    methods
       function obj = RayTracing(msh)
          obj.surf = msh; % should be a 2D mesh
          nNodes = size(msh.surfaces,2);
          surfCentroid = zeros(msh.nSurfaces,3);
          % compute surface centroids
          for i = 1:3
             tmp = msh.coordinates(msh.surfaces',i);
             tmp = reshape(tmp,nNodes,[]);
             surfCentroid(:,i) = sum(tmp,1)/nNodes;
          end
          obj.surfCentroids = surfCentroid;
            [obj.BBtree,obj.treeNodes, obj.leaf2elem] = obj.buildTree(msh);
        end

        function [BBtree, treeNodes, leaf2elem] = buildTree(obj,msh)
            % build the hierarchical tree of bounding boxes
            BBtree = zeros(2*msh.nSurfaces-1, 2);
            treeNodes = zeros(2*msh.nSurfaces-1,6);
            % store root polytop
            elemMap = false(msh.nSurfaces, size(BBtree,1));
            elemMap(:,1) = 1;
            leaf2elem = zeros(size(BBtree,1),1);
            k = 1;
            for i = 1:size(BBtree,1)
                % get elements of i-th node of the bounding volume tree
                surfList = find(elemMap(:,i));
                [treeNodes(i,:), leftCells, rightCells] = populateTree(obj,msh,surfList);
                if ~isempty(leftCells) % TreeNode is not a leaf Node
                    BBtree(i,:) = [2*k 2*k+1];
                    elemMap(leftCells,2*k) = 1;
                    elemMap(rightCells,2*k+1) = 1;
                    k = k + 1;
                else % Tree node is a leaf
                    leaf2elem(i) = surfList;
                end
            end
            obj.leaf2elem = leaf2elem;
        end

        function d = rayTrace(obj)
            % Build connectivity matrix between elements of 2 non
            % conforming grids
            % Uses Tandem traversal algorithm, starting from root of the
            % two trees
            d = tandemTraversal(obj,1,1);
        end

        function d = tandemTraversal(obj,t,origin,dir)
            % efficient algorithm to find element connectivity between the
            % two grids
            if ~rayTraceBox(obj,t,origin,dir)
                % Quit the procedure if the two bounding volumes do not
                % intersect
                return
            end

            if all(obj.BBtree(t,:) == 0) 
               % perform ray tracing on leaf cell
                d = intersectQuad(obj,t,origin,dir);
                return
            end

            if ~all(obj.BBtree(t,:) == 0)
                obj.tandemTraversal(obj.BBtree(t,1),origin,dir);
                obj.tandemTraversal(obj.BBtree(t,2),origin,dir);
            end
        end

        function [ktopVals, lCells, rCells] = populateTree(obj, msh, surfID)
            % INPUT: set of cells indices
            % OUTPUT: primitives of the bounding box
            %         Cells belonging to left and right child (if any)
            % get unique set of nodes belonging to input cells
            nodes = unique(msh.surfaces(surfID,:));
            % Store coordinates depending on 2D or 3D cases
            coords = msh.coordinates(nodes,1:3);
            ktopVals = [min(coords) max(coords)];
                        
            if length(surfID) > 1
                % split using cutting plane 
                % oriented like axis i 
                % passing trough point m
                [~,i] = max(abs(ktopVals(1:3) - ktopVals(4:6)));
                m = median(coords(:,i));
                surfPrim = obj.surfCentroids(surfID,i);
                id = surfPrim < m;
                lCells = surfID(id);
                rCells = surfID(~id);
                assert(length(lCells)+length(rCells) == length(surfID), 'Some elements left out from splitting procedure');
            else
                ktopVals = (ktopVals(:))';
                lCells = [];
                rCells = [];
            end
        end

        function id = rayTraceBox(obj,t,o,d)
           id = false;
           box = obj.BBtree(t,:); %[xmin xmax ymin ymax zmin zmax]
           tmin = zeros(3,1);
           tmax = zeros(3,1);
           for i = 1:3
              t1 = (box(2*i-1)-o(i))/d(i);
              t2 = (box(2*i)-o(i))/d(i);
              t = sort([t1 t2]);
              if d(i)~=0
                 tmin(i) = t(1); tmax(i) = t(2);
              else
                 tmin(i) = -Inf; tmax(i) = Inf;
              end
           end
           tmin = max(tmin); tmax = min(tmax);
           if tmin<=tmax && tmax>0
              id = true;
           end
        end

        function d = intersectQuad(obj,t,o,d)
           % get distance between a point and a mesh along a search
           % direction

           % d = nan : no intersection found with mesh
           % get element corresponding to leaf t
           el = obj.leaf2elem(t);
           % split quadrilateral in 2 triangles
           tri1 = obj.surf.coordinates(obj.surf.surfaces(el,1:3),:);
           tri2 = obj.surf.coordinates(obj.surf.surfaces(el,[1 3 4]),:);
           % perform ray tracing on each triangular subdivision
           d1 = MollerTrumbore(obj,o,d,tri1);
           d2 = MollerTrumbore(obj,o,d,tri2);
           d = min([d1 d2]);
        end

        function d = MollerTrumbore(obj,o,d,tri)
           [v1,v2,v3] = deal(tri(1,:)',tri(2,:)',tri(3,:)') ; 
           n = cross(v2-v1,v3-v1);
           if d'*n < 1e-4
              d = Inf;
              return
           else
              x = [d v2-v1 v3-v1]\(o-v1);
              if ~all([x(2)>=0 x(2)<=1 x(3)>=0 x(3)<=1])
                 d = Inf;
              else
                 d = x(1);
              end
           end
        end

    end
end


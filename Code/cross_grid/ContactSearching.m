classdef ContactSearching < handle
    %CONTACTSEARCHING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        msh1
        msh2
        BBtree1        % binary tree storying different mesh subidivisions
        BBtree2 
        polytop        % k-top used to bound a finite element
        treeNodes1     % For each node of the tree, store primitives along k/2 axis
        treeNodes2     % For each node of the tree, store primitives along k/2 axis
        leaf2elem1
        leaf2elem2
        elemConnectivity
        

    end
    
    methods
        function obj = ContactSearching(msh1,msh2,k)
            %CONTACTSEARCHING Construct an instance of this class
            %   Detailed explanation goes here
            obj.msh1 = msh1;
            obj.msh2 = msh2;
            obj.elemConnectivity = sparse(obj.msh1.nSurfaces, obj.msh2.nSurfaces);
            switch k
                case 8
                    obj.polytop = [1 0 -1 1;
                                   0 1 1 1];
                otherwise
                    error('%i - top discrete polytop is not supported',k);
            end
            [obj.BBtree1,obj.treeNodes1, obj.leaf2elem1] = obj.buildBBtree(msh1);
            [obj.BBtree2,obj.treeNodes2, obj.leaf2elem2] = obj.buildBBtree(msh2);
            obj.contactSearch();
        end

        function [BBtree, treeNodes, leaf2elem] = buildBBtree(obj,msh)
            % build the hierarchical tree of bounding boxes
            BBtree = zeros(2*msh.nSurfaces-1, 2);
            treeNodes = zeros(2*msh.nSurfaces-1, 2*size(obj.polytop,2));
            % store root polytop
            elemMap = sparse(msh.nSurfaces, size(BBtree,1));
            elemMap(:,1) = 1;
            leaf2elem = zeros(size(BBtree,1),1);
            k = 1;
            for i = 1:size(BBtree,1)
                % get elements of i-th node of the bounding volume tree
                surfList = find(elemMap(:,i));
                [treeNodes(i,:), leftCells, rightCells] = populateTreeNode(obj,msh,surfList);
                if ~isempty(leftCells) % TreeNode is not a leaf Node
                    BBtree(i,:) = [2*k 2*k+1];
                    elemMap(leftCells,2*k) = 1;
                    elemMap(rightCells,2*k+1) = 1;
                    k = k + 1;
                else % Tree node is a leaf
                    leaf2elem(i) = surfList;
                end
            end          
        end

        function contactSearch(obj)
            % Build connectivity matrix between elements of 2 non
            % conforming grids
            % Uses Tandem traversal algorithm, starting from root of the
            % two trees
            tandemTraversal(obj,1,1)
        end

        function tandemTraversal(obj, tNode1, tNode2)
            % efficient algorithm to find element connectivity between the
            % two grids
            if ~checkIntersection(obj, tNode1, tNode2)
                % Quit the procedure if the two bounding volumes do not
                % intersect
                return
            end

            if all(obj.BBtree1(tNode1,:) == 0) && all(obj.BBtree2(tNode2,:) == 0)
                obj.elemConnectivity(obj.leaf2elem1(tNode1),obj.leaf2elem2(tNode2)) = 1;
                return
            end

            if ~all(obj.BBtree1(tNode1,:) == 0)
                obj.tandemTraversal(obj.BBtree1(tNode1,1), tNode2);
                obj.tandemTraversal(obj.BBtree1(tNode1,2), tNode2);
            else
                obj.tandemTraversal(tNode1, obj.BBtree2(tNode2,1));
                obj.tandemTraversal(tNode1, obj.BBtree2(tNode2,2));
            end
        end

        function [ktopVals, lCells, rCells] = populateTreeNode(obj, msh, surfID)
            % INPUT: set of cells indices
            % OUTPUT: k primitives defining the bounding polytop of the
            %         Cells belonging to left and right child (if any)

            % get unique set of nodes belonging to input cells
            nodes = unique(msh.surfaces(surfID,:));
            if size(obj.polytop,2) == 4
                % assuming 2D case, discard z coordinate
                coords = msh.coordinates(nodes,1:2);
            end
            prim = coords*obj.polytop;
            ktopVals = [min(prim); max(prim)];
            % reduce slightly the size of the k-top (easy in 2D)
            red = min(abs(ktopVals(1,:) - ktopVals(2,:)));
            ktopVals = ktopVals + [0.01*red; -0.01*red];
            
            
            if length(surfID) > 1
                % split using cutting plane 
                % oriented like axis i 
                % passing trough point m
                [~,i] = max(abs(ktopVals(1,:) - ktopVals(2,:)));
                m = median(prim(:,i));
                ktopVals = (ktopVals(:))';
                surfPrim = msh.surfaceCentroid(surfID,1:2)*obj.polytop(:,i);
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

        function out = checkIntersection(obj, t1, t2)
            % check intersection between 2 polytops
            % if any of the primitives interval do not intersect, then
            % there's no intersection
            ktop1 = obj.treeNodes1(t1,:);
            ktop2 = obj.treeNodes2(t2,:);
            ktop1 = reshape(ktop1,2,[]);
            ktop2 = reshape(ktop2,2,[]);
            out = ~any([ktop1(1,:) > ktop2(2,:), ktop2(1,:) > ktop1(2,:)]);            
        end
    end
end


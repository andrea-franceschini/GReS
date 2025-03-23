classdef postProc
    % General purpose post processing utilities
    % Error analysis
    % Paraview plotting
    
    properties
        mesh
        element
        u
        u_ex
    end
    
    methods
        function obj = postProc(msh, u, u_ex, varargin)
            obj.mesh = msh; 
            if isempty(varargin)
                obj.element = Elements(msh);
            else
                obj.element = Elements(msh,varargin{1});
            end 

            % varargin: subsets of nodes for error computation
            assert(length(u)==obj.mesh.nNodes, ['Approximate Solution array size' ...
                'not matching number of mesh nodes']);
            assert(length(u_ex)==obj.mesh.nNodes, ['Approximate Solution array size' ...
                'not matching number of mesh nodes']);

            obj.u = u;
            obj.u_ex = u_ex;
        end
        
        function L2 = computeL2error(obj, varargin)

            if ~isempty(varargin)
                nodesList = varargin{1};
            else
                nodesList = 1:obj.mesh.nNodes;
            end
            L2 = (obj.u_ex(nodesList) - obj.u(nodesList)).^2;
            L2(isinf(L2)) = 0;
            L2(isnan(L2)) = 0;
            nodM = computeNodeMeasure(obj);
            L2 = sqrt(L2'*nodM);
        end

        function semiH1 = computeSemiH1error(obj, varargin)
            % varargin: subsets of nodes for error computation
            if ~isempty(varargin)
                elemList = varargin{1};
            else
                if (isempty(obj.mesh.nCells)) || (obj.mesh.nCells == 0)
                    elemList = 1:obj.mesh.nCells;
                else
                    elemList = 1:obj.mesh.nSurfaces;
                end
            end

            semiH1 = zeros(length(elemList),1);
            for el = elemList
                if (isempty(obj.mesh.nCells)) || (obj.mesh.nCells == 0)
                    switch obj.mesh.surfaceVTKType(el)
                        case 5 % Triangle
                            top = obj.mesh.surfaces(el, 1:obj.mesh.surfaceNumVerts(el));
                            N = getDerBasisF(obj.element.tri,el);
                            [vol,~] = findAreaAndCentroid(obj.element.tri,el);
                            semiH1(el) = (N*(obj.u_ex(top)-obj.u(top)))'*(N*(obj.u_ex(top)-obj.u(top)));
                            semiH1(el) = semiH1(el)*vol;
                        case 9 % Quadrilateral
                            top = obj.mesh.surfaces(el, :);
                            [N,dJWeighed] = obj.element.quad.getDerBasisFAndDet(el,1);
                            Nu_trans = pagemtimes(N,(obj.u_ex(top)-obj.u(top)));
                            Hs = pagemtimes(Nu_trans,'ctranspose',Nu_trans,'none');
                            Hs= Hs.*reshape(dJWeighed,1,1,[]);
                            semiH1(el) = sum(Hs,3);
                    end
                else
                    switch obj.mesh.cellVTKType(el)
                        case 10 % Tetrahedra
                            top = masterMesh.cells(el, 1:obj.mesh.cellNumVerts(el));
                            N = getDerBasisF(obj.element.tetra,el);
                            vol = findVolume(obj.element.tetra,el);
                            semiH1(el) = (N*(obj.u_ex(top)-obj.u(top)))'*(N*(obj.u_ex(top)-obj.u(top)));
                            semiH1(el) = semiH1(el)*vol;
                        case 12 % Hexahedra
                            top = obj.mesh.cells(el, :);
                            [N,dJWeighed] = obj.element.hexa.getDerBasisFAndDet(el,1);
                            Nu_trans = pagemtimes(N,(obj.u_ex(top)-obj.u(top)));
                            Hs = pagemtimes(Nu_trans,'ctranspose',Nu_trans,'none');
                            Hs= Hs.*reshape(dJWeighed,1,1,[]);
                            semiH1(el) = sum(Hs,3);
                    end
                end
            end
            semiH1 = sqrt(sum(semiH1));
        end



        function H1 = computeH1error(obj, varargin)
            % varargin: subsets of cells in the domain

            % define nodes and element subset
            if ~isempty(varargin)
                elemList = varargin{1};
                if (isempty(obj.mesh.nCells)) || (obj.mesh.nCells == 0)
                    nodesList = unique(obj.mesh.surfaces(elemList,:));
                else
                    nodesList = unique(obj.mesh.cells(elemList,:));
                end
            else
                nodesList = 1:obj.mesh.nNodes;
                if isempty(obj.mesh.nCells) || (obj.mesh.nCells == 0)
                    elemList = 1:obj.mesh.nSurfaces;
                else
                    elemList = 1:obj.mesh.nCells;
                end
            end
            % compute L2 and H1 seminorm errors
            L2 = computeL2error(obj, nodesList);
            semiH1 = computeSemiH1error(obj, elemList);
            H1 = sqrt(L2^2+semiH1^2);
        end

        function err_rel = computeRelError(obj)
            err_rel = abs((obj.u_ex- obj.u)./obj.u_ex);
            err_rel(isinf(err_rel)) = 0;
            err_rel(isnan(err_rel)) = 0;
            err_rel(err_rel > 1e4) = 0;
        end

        function nodMeasure = computeNodeMeasure(obj)
            % general class to compute the measure associated to each
            % node in 2D or 3D grids
            if (isempty(obj.mesh.nCells)) || (obj.mesh.nCells == 0)
                nodMeasure = obj.computeAreaNod();
            else
                nodMeasure = obj.computeVolNod();
            end
        end

        function gridSize = getGridSize(obj)
            gridSize = max(computeNodeMeasure(obj));
        end

        function a = computeAreaNod(obj)
            a = zeros(obj.mesh.nNodes,1);
            for el  = 1:obj.mesh.nSurfaces
                nodes = obj.mesh.surfaces(el,:);
                switch obj.mesh.surfaceVTKType(el)
                    case 5
                       [A,~] = obj.element.tri.findAreaAndCentroid(el);
                        a(nodes) = a(nodes) + A/3;
                    case 9
                        a(nodes) = a(nodes) + obj.element.quad.findNodeArea(el);
                end
            end
        end

        function vol = computeVolNod(obj)
            vol = zeros(obj.mesh.nNodes,1);
            for el  = 1:obj.mesh.nCells
                nodes = obj.mesh.cells(el,:);
                switch obj.mesh.cellVTKType(el)
                    case 10
                        vol(nodes) = vol(nodes) + obj.element.tetra.findVolume(el)/4;
                    case 12
                        vol(nodes) = vol(nodes) + obj.element.hexa.findNodeVolume(el);
                end
            end
        end
    end
end


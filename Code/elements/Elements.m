classdef Elements < handle
    % ELEMENT General element class

    properties (Access = public)
        % Number of elements by type
        nCellsByType   % nCellsByType = [#tetra, #hexa, #wed, #pyr]
        nSurfByType    % nCellsByType = [#tetra, #hexa, #wed, #pyr]
        nNodesElem = [4, 8, 6, 5]
        cellCentroid
        vol
        tetra
        hexa
        tri
        quad
        quadL
        indB
        indB2D
        mesh
    end

    properties (Access = private)
        % Element type
        %     elemType;
        % Mesh object
        %mesh
        GaussPts = []
    end

    methods (Access = public)
        % Class constructor method
        function obj = Elements(mesh,varargin)
            % Calling the function to set element data
            data = varargin;
            nIn = nargin;
            obj.setElementData(nIn,mesh,data);
            obj.computeCellProperties;
        end

        function computeCellProperties(obj)
            idCell = 1:obj.mesh.nCells;
            idTetra = idCell(obj.mesh.cellVTKType == 10);
            idHexa = idCell(obj.mesh.cellVTKType == 12);
            if ~isempty(idTetra)
                obj.vol(idTetra) = obj.tetra.findVolume(idTetra);
                obj.cellCentroid(idTetra,:) = computeCentroidGeneral(obj,idTetra);
            end
            if ~isempty(idHexa)
                [obj.vol(idHexa),obj.cellCentroid(idHexa,:)]= obj.hexa.findVolumeAndCentroid(idHexa);
                %          obj.vol(idHexa) = tmpVol;
                %          obj.cellCentroid(idHexa,:) = tmpCentroid;
            end
        end
    end

    methods (Access = private)
        % Assigns element data (properties of the 'Mesh' object)
        function setElementData(obj,nIn,mesh,data)
            obj.mesh = mesh;
            if nIn > 1
                obj.GaussPts = data{1};
            end
            %
            obj.nCellsByType = histc(obj.mesh.cellVTKType,[10, 12, 13, 14]);
            % in some cases histc may produce an empty array
            if isempty(obj.nCellsByType)
                obj.nCellsByType = zeros(4,1);
            end
            %
            if obj.nCellsByType(1) > 0
                obj.tetra = Tetrahedron(obj.mesh);
            end
            if obj.nCellsByType(2) > 0
                obj.hexa = Hexahedron(obj.mesh,obj.GaussPts);
            end
            obj.vol = zeros(obj.mesh.nCells,1);
            obj.cellCentroid = zeros(obj.mesh.nCells,3);
            %
            if obj.nCellsByType(2) == 0
                l1 = 4;
            else
                l1 = 8*obj.GaussPts.nNode;
            end
            obj.indB = zeros(9*l1,2);
            obj.indB(:,1) = repmat([1, 2, 3, 2, 1, 3, 3, 2, 1],[1,l1]);
            obj.indB(:,2) = repmat([1, 4, 6, 8,10,11,15,17,18],[1,l1]);
            obj.indB(:,1) = obj.indB(:,1) + repelem(3*(0:(l1-1))',9);
            obj.indB(:,2) = obj.indB(:,2) + repelem(18*(0:(l1-1))',9);

            % 2D elements (triangles, quadrilateral)
            obj.nSurfByType = histc(obj.mesh.surfaceVTKType,[5, 9]);
            if obj.nSurfByType(1) > 0
                obj.tri = Triangle(obj.mesh);
            end
            if obj.nSurfByType(2) > 0
                if any(obj.mesh.surfaceNumVerts==4)
                    obj.quad = Quadrilateral(obj.mesh,obj.GaussPts);
                elseif any(obj.mesh.surfaceNumVerts==8)
                    obj.quad = Quad8(obj.mesh,obj.GaussPts);
                    if obj.mesh.cartGrid
                        obj.quadL = Quadrilateral(getQuad4mesh(obj.mesh),obj.GaussPts);
                    end
                end
            end

            if obj.nSurfByType(2) == 0
                l1 = 3;
            else
                l1 = 4*obj.GaussPts.nNode;
            end
            % Matrix of derivatives of the basis functions for 2D elements
            obj.indB2D = zeros(4*l1,2);
            obj.indB2D(:,1) = repmat([1, 2, 2, 1],[1,l1]);
            obj.indB2D(:,2) = repmat([1, 3, 5, 6],[1,l1]);
            obj.indB2D(:,1) = obj.indB2D(:,1) + repelem(2*(0:(l1-1))',4);
            obj.indB2D(:,2) = obj.indB2D(:,2) + repelem(6*(0:(l1-1))',4);
        end

        function tetraCtr = computeCentroidGeneral(obj,idTetra)
            tetraCtr = sparse(repelem(1:length(idTetra),obj.mesh.cellNumVerts(idTetra)), ...
                nonzeros((obj.mesh.cells(idTetra,:))'), ...
                repelem((obj.mesh.cellNumVerts(idTetra)).^(-1), ...
                obj.mesh.cellNumVerts(idTetra)),length(idTetra),obj.mesh.nNodes) ...
                * obj.mesh.coordinates;
            %     end
        end
        %     function computeTetraCentroid(obj)
        %       obj.cellCentroid = sparse(repelem(1:obj.mesh.nCells,obj.mesh.cellNumVerts), ...
        %           nonzeros((obj.mesh.cells)'),repelem((obj.mesh.cellNumVerts).^(-1),obj.mesh.cellNumVerts),obj.mesh.nCells,obj.mesh.nNodes) ...
        %           * obj.mesh.coordinates;
        %     end
    end

end

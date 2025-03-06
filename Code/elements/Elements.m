classdef Elements < handle
    % ELEMENT General element class

    properties (Access = public)
        % Number of elements by type
        nCellsByType   % nCellsByType = [#tetra, #hexa, #wed, #pyr]
        nSurfByType    % nCellsByType = [#tetra, #hexa, #wed, #pyr]
        tetra
        hexa
        tri
        quad
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
           % compute volume/area and centroids of cell/surfaces in the mesh
           idCell = 1:obj.mesh.nCells;
           idSurf = 1:obj.mesh.nSurfaces;
           idTetra = idCell(obj.mesh.cellVTKType == 10);
           idHexa = idCell(obj.mesh.cellVTKType == 12);
           idTri = idSurf(obj.mesh.surfaceVTKType == 5);
           idQuad = idSurf(obj.mesh.surfaceVTKType == 9);
           if ~isempty(idTetra)
              [obj.mesh.cellVolume(idTetra), obj.mesh.cellCentroid(idTetra,:)]  = ...
                 obj.tetra.findVolumeAndCentroid(idTetra);
           end
           if ~isempty(idHexa)
              [obj.mesh.cellVolume(idHexa), obj.mesh.cellCentroid(idHexa,:)]  = ...
                 obj.hexa.findVolumeAndCentroid(idHexa);
           end
           if ~isempty(idTri)
              [obj.mesh.surfArea(idTri), obj.mesh.surfCentroid(idTri,:)]  = ...
                 obj.tri.findAreaAndCentroid(idTri);
           end
           if ~isempty(idQuad)
              [obj.mesh.surfArea(idQuad), obj.mesh.surfCentroid(idQuad,:)]  = ...
                 obj.quad.findAreaAndCentroid(idQuad);
           end
        end
    end

    methods (Access = private)
        function setElementData(obj,nIn,mesh,data)
           % Create instances of element type classes
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
            %
            % 2D elements (triangles, quadrilateral)
            obj.nSurfByType = histc(obj.mesh.surfaceVTKType,[5, 9]);
            if obj.nSurfByType(1) > 0
                obj.tri = Triangle(obj.mesh,obj.GaussPts);
            end
            if obj.nSurfByType(2) > 0
                if any(obj.mesh.surfaceNumVerts==4)
                    obj.quad = Quadrilateral(obj.mesh,obj.GaussPts);
                elseif any(obj.mesh.surfaceNumVerts==8)
                   error('Only 4 nodes quadrilateral are implemented')
                end
            end
        end
    end

end
classdef Elements < handle
  % ELEMENTS General element manager

  properties (Access = public)
    % Mesh reference
    mesh

    % list of elements istances 
    % [#tri, #quad, #tetra, #hexa #wed, #pyr, #quad9, #hexa27]
    mapVTK2elem
    elems = cell(1,8)

    % Element and surface type counts
    nCellsByType = zeros(1,5)     % [#tetra, #hexa, #wed, #pyr]
    nSurfByType = zeros(1,3)      % [#tri, #quad]
  end

  properties (Constant)
    % available finite elements 
    vtk3DTypeList = [10,12,13,14,29]
    vtk2DTypeList = [5,9,28]
    nNodesElem = [3,4,9,...       % 2D elements
                  4,8,27]         % 3D elements
  end


  properties (Access = protected)
    Jref      % basis function ref. gradient
    Nref      % basis function ref.
    Jb        % bubble basis function ref. gradient
    Nb        % bubble basis function ref.
  end

  methods (Access = public)
    function obj = Elements(mesh,nGaussPt)
        obj.mesh = mesh;
        if nargin < 2
          nGaussPt = 1;
        end
        setElementData(obj,mesh,nGaussPt);
    end


    function createElement(obj,id,msh,ng)
      k = obj.mapVTK2elem(id);
      switch id
        case 12
          obj.elems{k} = Hexahedron(ng,msh);
        case 10
          obj.elems{k} = Tetrahedron(ng,msh);
        case 9
          obj.elems{k} = Quadrilateral(ng,msh);
        case 5
          obj.elems{k} = Triangle(ng,msh);
        case 28
          obj.elems{k} = QuadrilateralQuadratic(ng,msh);
        case 29
          obj.elems{k} = HexahedronQuadratic(ng,msh);
        otherwise
          error('Finite element of vtk type %i not yet implemented \n',id)
      end
      obj.elems{k}.computeProperties();
    end
  end

  methods (Access = private)
    function setElementData(obj,msh,g)

      obj.mapVTK2elem = zeros(max(obj.vtk3DTypeList),1);
      vtkList = [obj.vtk2DTypeList obj.vtk3DTypeList];
      obj.mapVTK2elem(vtkList) = 1:numel(obj.elems);

      % Create 3D finite elements in the mesh
      c3D = 0;
      for vtk = obj.vtk3DTypeList
        c3D = c3D+1;
        nC = sum(obj.mesh.cellVTKType == vtk);
        if nC > 0
          obj.nCellsByType(c3D) = nC;
          obj.createElement(vtk,msh,g);
        end
      end

      c2D = 0;
      % Create 2D finite elements in the mesh
      for vtk = obj.vtk2DTypeList
        c2D = c2D+1;
        nC = sum(obj.mesh.surfaceVTKType == vtk);
        if nC > 0
          obj.nSurfByType(c2D) = nC;
          obj.createElement(vtk,msh,g);
        end
      end
    end
  end


      % Count 2D surface types
%       obj.nSurfByType = histc(obj.mesh.surfaceVTKType, [5, 9]);
% 
%       if obj.nSurfByType(1) > 0
%         obj.tri = Triangle(obj.mesh, obj.GaussPts);
%       end
% 
%       if obj.nSurfByType(2) > 0
%         if any(obj.mesh.surfaceNumVerts == 4)
%           obj.quad = Quadrilateral(obj.mesh, obj.GaussPts);
%         elseif any(obj.mesh.surfaceNumVerts == 8)
%           obj.quad = Quad8(obj.mesh, obj.GaussPts);
%           if obj.mesh.cartGrid
%             obj.quadL = Quadrilateral(getQuad4mesh(obj.mesh), obj.GaussPts);
%           end
%         end
%       end

      % Build indB2D
%       if obj.nSurfByType(2) == 0
%         l1 = 3;
%       else
%         l1 = 4 * obj.GaussPts.nNode;
%       end
% 
%       obj.indB2D = zeros(4 * l1, 2);
%       obj.indB2D(:,1) = repmat([1, 2, 2, 1], [1, l1]);
%       obj.indB2D(:,2) = repmat([1, 3, 5, 6], [1, l1]);
%       obj.indB2D(:,1) = obj.indB2D(:,1) + repelem(2*(0:(l1-1))', 4);
%       obj.indB2D(:,2) = obj.indB2D(:,2) + repelem(6*(0:(l1-1))', 4);

  methods (Access = public)
    % these methods must are mandatory in each finite element class
    function element = getElement(obj,vtkId)
      element = obj.elems{obj.mapVTK2elem(vtkId)};
    end

    function element = getElementByID(obj,id)
      vtkId = obj.mesh.cellVTKType(id);
      element = obj.elems{obj.mapVTK2elem(vtkId)};
    end

    function areaNod = findNodeArea(obj,el)
      i = obj.mesh.surfaceVTKType(el);
      areaNod = findNodeArea(getElement(obj,i),el);
    end

    function volNod = findNodeVolume(obj,el)
      i = obj.mesh.cellVTKType(el);
      volNod = findNodeVolume(getElement(obj,i),el);
      volNod = reshape(volNod,[],1);
    end

    function N = getNumbCellData(obj,varargin)
      % return total number of gauss point in mesh
      % useful to initialize stress/strain variables
      N = 0;
      k = 0;
      if ~isempty(varargin)
        assert(nargin==2,'Wrong number of input arguments');
        elemList = varargin{1};
      end
      for vtkID = obj.vtk3DTypeList
        k = k+1;
        elem = getElement(obj,vtkID);
        if ~isempty(elem)
          ng = elem.GaussPts.nNode;
          if isempty(varargin)
            nc = obj.nCellsByType(k);
          else
            nc = sum(obj.mesh.cellVTKType(elemList) == vtkID);
          end
          N = N + ng*nc;
        end
      end
    end

%     function N = getNumbNodeData(obj,varargin)
%       % return number of nodes for each cell in input
%       N = 0;
%       k = 0;
%       if ~isempty(varargin)
%         assert(nargin==2,'Wrong number of input arguments');
%         elemList = varargin{1};
%       end
%       for vtkID = obj.vtk3DTypeList
%         k = k+1;
%         elem = getElement(obj,vtkID);
%         if ~isempty(elem)
%           if isempty(varargin)
%             nc = obj.nCellsByType(k);
%           else
%             nc = sum(obj.mesh.cellVTKType(elemList) == vtkID);
%           end
%           N = N + elem.nNode*nc;
%         end
%       end
%     end
  end

end

classdef (Abstract) MortarQuadrature < handle
  % Mortar integration on [interface, slave, master] neighbor triplets.

  properties (SetAccess = private)
    grid
    gpCoords = cell(0,3)   % reference coordinates, one row per triplet
    weights = cell(0,1)    % physical integration weights per triplet
    multiplierType
    gaussOrder
  end

  properties (SetAccess = private)
    areaSlave = []
  end


  properties (Access = private)
    % store in memory transform operators for dual basis construction
    dualTransforms
  end

  methods (Abstract)

    [xi,w] = processMortarPair(obj,tripletId,elements)


  end

  methods

    function obj = MortarQuadrature(multType,interfGrid,input)

      options = readInput(struct('gaussOrder',4),input);
      obj.gaussOrder = options.gaussOrder;

      obj.multiplierType = validatestring(multType,{'P0','standard','dual'});

      obj.grid = interfGrid;

    end


    function processMortarTriplets(obj)

      if ~isempty(obj.weights)
        return
      end

      % Segment preparation can prune and renumber the intersection triplets.
      obj.prepareQuadrature();
      elements = obj.createElements();

      neighbors = obj.grid.neighbors;
      n = size(neighbors,1);
      coords = cell(n,3);
      w = cell(n,1);
      for k = 1:n

        ids = neighbors(k,:);
        elems = {elements{1}{ids(1)},elements{2}{ids(2)},elements{3}{ids(3)}};

        [coords(k,:),w{k}] = obj.processMortarPair(k,elems);

      end

      % keep only active triplets (at least a gauss weight)
      active = ~cellfun(@isempty,w);
      obj.grid.removeNeighbors(~active);
      obj.dualTransforms = cell(obj.grid.cells.num,1);
      obj.gpCoords = coords(active,:);
      obj.weights = w(active);
      obj.areaSlave = accumarray(obj.grid.neighbors(:,2),cellfun(@sum,obj.weights), ...
        [obj.grid.traceGrids(1).cells.num,1]);

    end

    function ids = getActiveTriplets(obj)
      ids = (1:numel(obj.weights)).';
    end

    function [xiInterf,xiSlave,xiMaster] = getGPCoords(obj,triplet)
      [xiInterf,xiSlave,xiMaster] = obj.gpCoords{triplet,:};
    end

    function xi = getSlaveGPCoords(obj,triplet)
      xi = obj.gpCoords{triplet,2};
    end

    function xi = getMasterGPCoords(obj,triplet)
      xi = obj.gpCoords{triplet,3};
    end

    function w = getIntegrationWeights(obj,triplet)
      w = obj.weights{triplet};
    end

    function [Nslave,Nmaster,Nmult] = getMortarBasisFunctions(obj,triplet,elemInterf,elemMaster,elemSlave)

      % return basis functions for a given triplet of connected elmements
      ids = obj.grid.neighbors(triplet,:);
      [xiInterf,xiSlave,xiMaster] = obj.getGPCoords(triplet);

      Nslave = elemSlave.computeBasisF(xiSlave);
      Nmaster = elemMaster.computeBasisF(xiMaster);

      Nmult = obj.computeMultiplierBasis(ids(1),elemInterf,xiInterf);

    end


    function Nmult = computeMultiplierBasis(obj,elId,elem,xi)

      if strcmp(obj.multiplierType,'P0')
        Nmult = ones(size(xi,1),1);
        return
      end

      % standard case
      Nmult = elem.computeBasisF(xi);

      if strcmp(obj.multiplierType,'dual')


      error('Temporarily unavailable')


        % store dual mapping
        if isempty(obj.dualTransforms{elId})
          % Biorthogonality is defined over the whole interface facet.
          % Use a rule that integrates products of the facet basis functions.
          rule = FiniteElementType.create(elem.vtkType,obj.grid, ...
            'gaussOrder',max([2,2*elem.minGaussOrder-1,obj.gaussOrder]));
          N = rule.getBasisFinGPoints();
          [~,w] = rule.getDerBasisFAndDet(elId);
          w = sign(sum(w))*w(:);
          M = N'*(N.*w);
          D = diag(N'*w);
          obj.dualTransforms{elId} = M\D;
        end

        Nmult = Nmult*obj.dualTransforms{elId};
      end
    end
  end


  methods (Access = protected)
    function prepareQuadrature(~)
    end
  end

  methods (Access = private)

    function elements = createElements(obj)

      % helper function allocating ready to use elements for mortar
      % processing
      grids = {obj.grid,obj.grid.traceGrids(1),obj.grid.traceGrids(2)};
      elements = cell(3,1);
      for side = 1:3
        mesh = grids{side};
        ids = double(mesh.cells.VTKType);
        types = cell(1,max([0;ids(:)]));
        for vtk = reshape(unique(ids),1,[])
          types{vtk} = FiniteElementType.create(vtk,mesh,'gaussOrder',obj.gaussOrder);
        end
        elements{side} = types(ids);
      end
    end
  end

  methods (Static, Access = protected)

    function basis = facetBasis(coords)
      basis = eye(2);
      if size(coords,2) == 3
        if size(coords,1) == 9
          coords = coords(1:4,:);
        end
        normal = computePolygonNormal(coords);
        [~,d1,d2] = pointToSurfaceProjection(coords(1,:),normal,coords);
        basis = [d1,d2];
      end
    end

    function [xi,converged] = mapToReference(elem,coords,points)
      % Invert a facet map in the two-dimensional projection plane.
      n = size(points,1);
      xi = repmat(elem.centroid,n,1);
      converged = false(n,1);
      tolerance = 1e-10*max(vecnorm(coords - coords(1,:),2,2));

      for g = 1:n
        for iter = 1:20
          residual = points(g,:) - elem.computeBasisF(xi(g,:))*coords;
          if norm(residual) <= tolerance
            converged(g) = true;
            break
          end
          J = elem.computeDerBasisF(xi(g,:))*coords;
          if rcond(J) <= eps
            break
          end
          % With row coordinates, the Newton increment satisfies dXi*J = residual.
          xi(g,:) = xi(g,:) + residual/J;
        end
      end
    end
  end

  methods (Static)
    function mat = integrate(func,varargin)
      dJw = varargin{end};
      mat = func(varargin{1:end-1});
      mat = mat.*reshape(dJw,1,1,[]);
      mat = sum(mat,3);
    end
  end
end

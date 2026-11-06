classdef PeacemanWells < handle

  % Peaceman well models

  properties (SetAccess = protected)

    wells     % struct of wells defined using MRST addWells
    grid
    mat 
  

  end

  methods (Access = public)

    
    function obj = PeacemanWells(domain)

      obj.grid = domain.grid;
      obj.mat = domain.materials;

    end


    function addWell(obj,cellList,varargin)

      wId = numel(obj.wells) + 1;

      numC = numel(cellList);

      % define default structure
      default = struct(                               ...
        'dir'         , 'z',                          ...
        'name'        , sprintf('W%d', wId + 1),      ...
        'radius'      , 0.1,                          ...
        'control'     , "bhp",                        ...
        'time'        , 0.0,                          ...
        'sign'        , [],                           ...
        'value'         , 0,                            ...
        'WI'          , repmat(-1, [numC, 1]),        ...
        'skin'        , zeros(numC, 1),               ...
        'comp'        , [0.5 ,0.5],                        ...
        'refDepth'    , max(obj.grid.coordinates(:,3)));

  
      newWell = readInput(default,varargin{:});

      newWell.cellList = cellList(:);

      t = newWell.time;
      v = newWell.value;
      d = find(strcmp(["x","y","z"], newWell.dir));

      if numel(t) ~= numel(v)
        error('Well time schedule instants must coincide with value')
      end

      if isscalar(t)
        newWell.time = [t t+1];
        newWell.value = [v v];
      end

      if ~isscalar(newWell.dir)
        error("Direction of the well must be a scalar")
      end

      newWell.dir = d;

      if newWell.sign ~= -1 && newWell.sign ~= +1
        error('Wrong sign parameter for well %s. Expected +1 for injector and -1 for producer',newWell.name)
      end


      obj.wells = [obj.wells; newWell];

    end



    function computeWellProperties(obj,perm)


      for i = 1:numel(obj.wells)

        cellId = obj.wells(i).cellList;

        % extract diagonal components of well cell permeabilities
        Ki = perm(cellId,[1 5 9]);

        Hi = getCellDimensions(obj,cellId);

        obj.wells(i).effRadius = obj.computeEffectiveRadius(obj.wells(i).dir,Ki,Hi);

        obj.wells(i).WI = obj.computeWellIndex(obj.wells(i),obj.wells(i).dir,Ki,Hi);

        obj.wells(i).dZ = obj.wells(i).refDepth - obj.grid.cells.center(cellId,3);


      end
    end

    function val = getValue(obj,wId,t)

      W = obj.wells(wId);
      T = W.time;
      V = W.value;

      t = max(T(1), min(t, T(end)));   % Clamp to the available interval

      i = find(T <= t, 1, 'last');
      i = min(i, numel(T)-1);
      j = i + 1;

      alpha = (t - T(i)) / (T(j) - T(i));
      alpha = max(0, min(alpha, 1));

      val = (1-alpha)*V(i) + alpha*V(j);

    end


  end




    






  methods (Access = private)

    function h = getCellDimensions(obj,cellList)

      % size of a bounding box passing trough the face centroid of the
      % cells
      g = obj.grid;
      faceTopol = g.cells.cells2faces.getRowsMatrix(cellList);
      C = g.faces.center;
      h = zeros(numel(cellList),3);

      for i = 1:3
        c = C(:,i);
        ci = c(faceTopol);
        ci = reshape(ci,numel(cellList),[]);
        h(:,i) = max(ci,[],2) - min(ci,[],2);
      end

    end



    % function cells = getPerforatedCells(obj,wId)
    % 
    %   if nargin == 1
    %     cells = cell2mat(obj.cellList(:));
    %   else
    %     cells = obj.cellList{wId}(:);
    %   end
    % 
    % end

    end



  methods (Static)

    function r = computeEffectiveRadius(d,K,h)
      % ref: White et al. 2019

      v = 1:3;
      dd = v(v ~= d);

      % get permeabiliy

      h1 = h(:,dd(1));
      h2 = h(:,dd(2));

      r1 = K(:,dd(1))./K(:,dd(2));
      r2 = 1./r1;

      r = 0.28 * ((r1.^(0.5) .* h2.^2 + r2.^(0.5) .* h1.^2).^0.5)./(r1.^0.25 + r2.^0.25);

    end


    function WI = computeWellIndex(w,d,K,H)
      % ref: White et al. 2019

      rw = w.radius;
      sk = w.skin;
      r = w.effRadius;

      v = 1:3;
      dd = v(v ~= d);

      kk = K(:,dd(1)).*K(:,dd(2));

      WI = 2*pi*H(:,d).*sqrt(kk)./(log(r./rw) + sk);


    end


  end



end

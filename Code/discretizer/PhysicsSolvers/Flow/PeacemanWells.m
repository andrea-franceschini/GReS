classdef PeacemanWells < handle

  % Peaceman well models

  % Currently implemented for TPFA SinglePhaseFlow discretization

  properties (SetAccess = protected)

    wells     % struct of wells defined using MRST addWells
    cellList
    domain 
    ADI
  

  end

  methods (Access = public)

    
    function obj = PeacemanWells(domain)

      obj.domain = domain;
      obj.cellList = cell([]);

    end


    function addWell(obj,cellList,varargin)

      wId = numel(obj.wells) + 1;

      numC = numel(cellList);

      obj.cellList{wId} = cellList(:);

      grid = obj.domain.grid;


      % define default structure
      default = struct(                               ...
        'dir'         , 'z',                          ...
        'name'        , sprintf('W%d', wId + 1),      ...
        'radius'      , 0.1,                          ...
        'control'     , "bhp",                        ...
        'time'        , 0.0,                          ...
        'val'         , 0,                            ...
        'WI'          , repmat(-1, [numC, 1]),        ...
        'skin'        , zeros(numC, 1),               ...
        'refDepth'    , max(grid.coordinates(:,3)));

  
      newWell = readInput(default,varargin{:});

      t = newWell.time;
      v = newWell.value;

      if numel(t) ~= numel(v)
        error('Well time schedule instants must coincide with value')
      end

      if isscalar(t)
        newWell.time = [t t+1];
        newWell.value = [v v];
      end

      d = find(strcmp(["x","y","z"], newWell.dir));

      if ~isscalar(d)
        error("Direction of the well must be a scalar")
      end

      % get permeability matrix
      K = zeros(numC,9);
      ctags = obj.domain.grid.cells.tag(cellList);

      mat = obj.domain.materials;

      for i = reshape(unique(ctags),1,[])
        id = ctags == i;
        k = mat.getMaterial(i).PorousRock.getPermVector();
        K(id,:) = repmat(k,sum(id),1);
      end

      % extract diagonal components
      K = K(:,[1 5 9]);

      H = getCellDimensions(obj,cellList);

      newWell.effRadius = obj.computeEffectiveRadius(d,K,H);

      newWell.WI = obj.computeWellIndex(newWell,d,K,H);

      newWell.dZ = grid.cells.center(:,3)

      obj.wells = [obj.wells; newWell];

    end




    function assembleWells(obj)

      % assemble well contribution within linear system



      
    end


    function initialize(obj)

      % setup the well ADI equations

      nw = numel(obj.wells);

      WI = vertcat(obj.wells(:).WI);    % well-indices
      dz = vertcat(obj.wells(:).dZ);    % connection depth relative to bottom-hole

      % using MRST ADI framework
      f = obj.domain.materials.getFluid;
      gamma = f.getSpecificWeight;
      mu = f.getDynamicViscosity;
      % applies only to perforated cells
      obj.ADI.fluxEq = @(p,bhp) WI .* (1 / mu) .* (bhp - gamma*dz - p(obj.getPerforatedCells));
      obj.ADI.p = initVariablesADI(getStateInit(obj).pressure);

      %rateEq = @(p,bhp,qS)  qS-sum(q_conn(p, bhp))/rhoS;
      %ctrlEq = @(bhp)       bhp-100*barsa;

      if any(strcmp([obj.wells.control],"rate"))

        error('Only BHP control is currently supported')

      end





        


      end

      % if flow control, we need one extra equation for each well,
      % corresponding to a new BHP variable that must be registered in the
      % model



    end







  methods (Access = private)


    function h = getCellDimensions(obj,cellList)

      % size of a bounding box passing trough the face centroid of the
      % cells
      g = obj.domain.grid;
      faceTopol = g.cells.cells2faces.getRowsMatrix(cellList);
      C = g.faces.center;
      h = zeros(numel(cellList),3);

      for i = 1:3
        c = C(:,i);
        ci = c(faceTopol);
        h(:,i) = max(ci,[],2) - min(ci,[],2);
      end

    end



    function cells = getPerforatedCells(obj,wId)

      if nargin == 1
        cells = cell2mat(obj.cellList(:));
      else
        cells = obj.cellList{wId}(:);
      end

    end

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

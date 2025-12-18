classdef SedimentsMap < handle
  % SEDIMENTSMAP Manages spatial/temporal sedimentation history and
  % interpolation.
  %
  %   Loads events from a file and calculates material accumulation over
  %   a time interval [t0, t0+dt]. Supports 'linear' (trapezoidal) and
  %   piecewise constant ('ramp') interpolation types.

  properties (Access = public)
    timeList = []
    dataStruct
  end

  properties (Access = private)
    type
    dim (1,2) uint64
    nmat uint16
  end

  methods
    function obj = SedimentsMap(data,nmat,dims)
      % Constructor: Initializes grid and loads event file.
      if nargin == 0
        return
      end
      obj.dim = dims;
      obj.nmat = nmat;
      if ~isfield(data,"file")
        error("Maps for the sedimentation missing!");
      end
      obj.constructor(data.file);
    end

    function map = getSedimentationMap(obj,t0,dt)
      % GETSEDIMENTATIONMAP Integrates sedimentation rate from t0 to t0+dt.
      %   Directs logic to specific interpolation methods based on 'obj.type'.
      % Returns the sediment map for time t0 and duration dt.

      ntimesteps = length(obj.timeList);
      pos = find(obj.timeList <= t0, 1, 'last');
      posF = find(t0+dt <= obj.timeList, 1, 'first');

      % Boundary checks
      if isempty(pos)
        warning("The sediment map was not defined! Returning zeros.");
        map = zeros(prod(obj.dim), obj.nmat);
        return;
      end
      if isempty(posF)
        posF=pos;        
      end

      % Routing to interpolation logic
      if ntimesteps > 1 && pos ~= posF
        if strcmpi(obj.type, "linear")
          map = obj.LinearInter(dt, pos, posF);
        else
          map = obj.RampInter(dt, pos, posF);
        end
      else
        % Single event or within the same time window
        map = obj.processEvent(pos) * dt;
      end
    end

  end

  methods (Access = private)
    function constructor(obj, file)
      % Parses file, sets type, and sorts events by time.
      listdata = readstruct(file, AttributeSuffix="");
      obj.type = lower(listdata.type);
      obj.dataStruct = listdata.Event;

      % Remove duplicate timestamps and sort chronologically
      [~, uniqueIdx] = unique([listdata.Event.time]);
      obj.dataStruct = listdata.Event(uniqueIdx);
      obj.timeList = [obj.dataStruct.time];
    end

    function map = processEvent(obj, pos)
      % Aggregates Uniform and Map data for a specific event index.
      map = zeros(prod(obj.dim), obj.nmat);
      checkMat = true(obj.nmat, 1);
      event = obj.dataStruct(pos);

      % Handle Uniform distributions
      if isfield(event, "Uniform")
        for item = event.Uniform
          if checkMat(item.materialFlag)
            map(:, item.materialFlag) = item.value;
            checkMat(item.materialFlag) = false;
          end
        end
      end

      % Handle File-based Maps
      if isfield(event, "Map")
        for item = event.Map
          if checkMat(item.materialFlag) && isfile(item.file)
            map(:, item.materialFlag) = obj.dataMap(item);
            checkMat(item.materialFlag) = false;
          end
        end
      end
    end

    function values = dataMap(obj, data)
      % Loads map from file and validates spatial dimensions.
      mapDim = sscanf(data.division, '%lu,%lu')';
      if isequal(obj.dim, mapDim)
        values = load(data.file);
      else
        % TODO : Create an interpolation for the map values
        warning("The sediment map has a mismatch dimension" + ...
          " with the domain. Returning zeros.");
        values = zeros(prod(obj.dim), 1);
      end
    end

    function map = LinearInter(obj, dt, pos, posF)
      % LINEARINTER Performs trapezoidal integration across event intervals.
      map = zeros(prod(obj.dim), obj.nmat);
      dtAcc = 0;
      for i=1:(posF-pos)
        tSt = obj.timeList(pos+i-1);
        tEd = obj.timeList(pos+i);
        v1 = obj.processEvent(pos+i-1);
        v2 = obj.processEvent(pos+i);

        if i==(posF-pos)
          dl = dt-dtAcc;
          % Interpolated rate at time t0
          % y = A x + B
          tDiff = 1/(tEd - tSt);
          Aterm = tDiff*(v2 - v1);
          Bterm = tDiff*(tEd*v1 - tSt*v2);
          v2 = Aterm*(tSt+dl)+Bterm;
        else
          dl = tEd-tSt;
          dtAcc = dtAcc + dl;
        end
        map = map + (dl/2)*(v1+v2);
      end
    end

    function map = RampInter(obj, dt, pos, posF)
      % RAMPINTER Performs piecewise constant integration (step function).
      map = zeros(prod(obj.dim), obj.nmat);
      dtAcc = 0;
      for i=1:(posF-pos)
        % Constant by part's between Event(pos) and Event(pos+1)
        tSt = obj.timeList(pos+i-1);
        tEd   = obj.timeList(pos+i);

        if i==(posF-pos)
          dl = dt-dtAcc;
        else
          dl = tEd-tSt;
          dtAcc = dtAcc + dl;
        end
        v = obj.processEvent(pos+i-1);
        map = map + v * dl;
      end
    end

  end
end
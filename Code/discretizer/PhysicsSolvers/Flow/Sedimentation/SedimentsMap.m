classdef SedimentsMap < handle
  % SEDIMENTATION MAP - General class to store the sedimentation history
  %   Detailed explanation goes here

  properties (Access = public)
    timeList = []
    map = []
    dataStruct
  end

  properties (GetAccess=public, SetAccess=private)
    type
  end

  properties (Access = private)
    dim (1,2) uint64
    nmat uint16
  end

  methods
    function obj = SedimentsMap(data,nmat,dims)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here

      % obj.db = containers.Map('KeyType','char','ValueType','any');
      if nargin == 0
        return
      end
      obj.dim = dims;
      obj.nmat = nmat;
      % obj.map = zeros(prod(dims),nmat);
      if ~isfield(data,"file")
        error("Maps for the sedimentation missing!");
      end
      obj.constructor(data.file);
    end

    function map = getSedimentationMap(obj,t0,dt)

      map = zeros(prod(obj.dim),obj.nmat);
      tf=t0+dt;
      
      ntimesteps = length(obj.timeList);
      posA=1;
      posB=1;
      if ntimesteps>1
        for t=1:ntimesteps
          tref=obj.timeList(t);
          if tref<t0
            posA=posA+1;
          end
          if tref<tf
            posB=posB+1;
          end

        end
      else
      end


    end


  end

  methods (Access = private)
    function constructor(obj,file)
      % Reading
      listdata=readstruct(file, AttributeSuffix="");
      if isfield(listdata,"type")
        obj.type = lower(listdata.type);
      else
        obj.type = "ramp";
      end
      obj.dataStruct = listdata.Event;

      % Constructing the time array.
      numEvent = length(listdata.Event);
      toDelete = false(numEvent,1);
      for t=1:numEvent
        tmp=listdata.Event(t);
        if ~ismember(tmp.time,obj.timeList)
          obj.timeList(end+1)=tmp.time;
        else
          toDelete(t)=true;
        end
      end

      % Deleting duplicates and Ordering as function of the time.
      obj.dataStruct(toDelete)=[];
      [~,idx] = sort(obj.timeList);
      obj.dataStruct=obj.dataStruct(idx);
      obj.timeList=obj.timeList(idx);
    end

    function processEvent()
    end

    function tmpEvent()

      data = listdata(1).Event;

      if ~isfield(data,"materialFlag")
        data = listdata(2);
      end

      if isfield(data,"uniform")
        obj.map(:,data.materialFlag) = data.uniform;
        return
      end

      if isfield(data,"division")
        tmp = str2num(data.division);
        % tmp
        obj.map(:,data.materialFlag) = data.uniform;
        return
      end


    end

    function constructorStruct_old(obj,data)
      % obj.dim = str2num(data.division);
      tmp = data.sediment_rate;
      tmp = str2num(tmp);
      obj.map = reshape(tmp,obj.dim);
    end

    

  end
end
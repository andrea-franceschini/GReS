classdef SedimentsMap < handle
  % SEDIMENTATION MAP - General class to store the sedimentation history
  %   Detailed explanation goes here

  properties (Access = private)
    % db 
    map = []
    dim (1,2) int64
    % t
  end

  methods
    function obj = SedimentsMap(varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here

      % obj.db = containers.Map('KeyType','char','ValueType','any');
      if nargin == 0
        return
      end
      obj.dim = varargin{1}(1:2);
      obj.constructorStruct(varargin{2});
    end

    function map = getSedimentationMap(obj,t)
      map = obj.map;
    end


  end

  methods (Access = private)
    function constructorStruct(obj,data)
      % obj.dim = str2num(data.division);
      tmp = data.sediment_rate;
      tmp = str2num(tmp);
      obj.map = reshape(tmp,obj.dim);
    end
  end
end
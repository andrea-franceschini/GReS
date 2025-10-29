classdef BoundaryEntities < handle
  % BOUNDARY CONDITIONS ENTITY class
  % Constrained entities IDs and values of the BCs are defined by the user.

  properties (Access = public)
    % Boundary condition identifier (mainly for error messages)
    name
    % Total number of constrained entities
    totEnts
    % Number of constrained entities for each degree of freedom
    nEntities
    % Indices of constrained degrees of freedom
    entities
    % Number of input times
    nTimes
    % Set of input times
    times
    % Set of files with boundary conditions in time
    dataFiles
    % Values of currently-stored boundary conditions
    availVals
    % Time id of currently-stored boundary conditions
    availSteps
  end
  
  properties (Access = private)
      dof
      load = false
  end

  methods (Access = public)
    % Class constructor method
    function obj = BoundaryEntities(name, setFile, times, dataFiles)
      % Calling the function to set object properties
      obj.setBC(name, setFile, times, dataFiles);
    end

    function vals = getValues(obj, t)
      % Load the boundary parameters
      [i1, i2] = bin_search(obj, t);
      if (~obj.load) & ~isequal(obj.availSteps,[i1;i2])
        obj.loadValues(t);
      end
      if (obj.nTimes == 1)
        vals = obj.availVals(:,1);
      else
        fac = (t - obj.times(1)) / (obj.times(2) - obj.times(1));
        vals = fac*(obj.availVals(:,2) - obj.availVals(:,1)) + obj.availVals(:,1);
      end
    end

    % function vals = getValues(obj, t)
    %   if (obj.nTimes == 1)
    %     vals = readDataSet(obj.dataFiles(1).fileName, obj.totEnts);
    %   else
    %     [i1, i2] = bin_search(obj, t);
    %     p1 = find(obj.availSteps == i1);
    %     p2 = find(obj.availSteps == i2);
    %     if (isempty(p1) && isempty(p2))
    %       p1 = 1;
    %       obj.availVals(:,1) = readDataSet(obj.dataFiles(i1).fileName, obj.totEnts);
    %       obj.availSteps(1) = i1;
    %       p2 = 2;
    %       obj.availVals(:,2) = readDataSet(obj.dataFiles(i2).fileName, obj.totEnts);
    %       obj.availSteps(2) = i2;
    %     elseif (~isempty(p1) && isempty(p2))
    %       p2 = 3 - p1;
    %       obj.availVals(:,p2) = readDataSet(obj.dataFiles(i2).fileName, obj.totEnts);
    %       obj.availSteps(p2) = i2;
    %     elseif (~isempty(p2) && isempty(p1))
    %       p1 = 3 - p2;
    %       obj.availVals(:,p1) = readDataSet(obj.dataFiles(i1).fileName, obj.totEnts);
    %       obj.availSteps(p1) = i1;
    %     end
    %     fac = (t - obj.times(i1)) / (obj.times(i2) - obj.times(i1));
    %     vals = fac*(obj.availVals(:,p2) - obj.availVals(:,p1)) + obj.availVals(:,p1);
    %   end
    % end

    function loadValues(obj, t)
      if (obj.nTimes == 1)
        obj.availVals(:,1) = readDataSet(obj.dataFiles(1).fileName, obj.totEnts);
      else
        [i1, i2] = bin_search(obj, t);
        p1 = find(obj.availSteps == i1);
        p2 = find(obj.availSteps == i2);
        if (isempty(p1) && isempty(p2))
          obj.availVals(:,1) = readDataSet(obj.dataFiles(i1).fileName, obj.totEnts);
          obj.availSteps(1) = i1;
          obj.availVals(:,2) = readDataSet(obj.dataFiles(i2).fileName, obj.totEnts);
          obj.availSteps(2) = i2;
        elseif (~isempty(p1) && isempty(p2))
          p2 = 3 - p1;
          obj.availVals(:,p2) = readDataSet(obj.dataFiles(i2).fileName, obj.totEnts);
          obj.availSteps(p2) = i2;
        elseif (~isempty(p2) && isempty(p1))
          p1 = 3 - p2;
          obj.availVals(:,p1) = readDataSet(obj.dataFiles(i1).fileName, obj.totEnts);
          obj.availSteps(p1) = i1;
        end
      end
      obj.load = true;
    end



  end

  methods (Access = private)
    function setBC(obj, name, setFile, times, dataFiles)
      obj.name = name;
      obj.times = times;
      obj.nTimes = length(times);
      obj.dataFiles = dataFiles;
      [obj.nEntities, obj.entities] = readEntitySet(setFile);
      obj.totEnts = sum(obj.nEntities);
      if (obj.totEnts == 0)
        error('No boundary conditions are prescribed for %s BC', name);
      end
      obj.availVals = zeros(obj.totEnts,2);
      obj.availSteps = zeros(2,1);
    end

    function [i1, i2] = bin_search(obj, t)
      i1 = 1;
      i2 = obj.nTimes;
      pos = floor(i2/2);
      while (i2-i1 > 1)
        if (obj.times(pos) >= t)
          % Left interval
          i2 = pos;
        else
          % Right interval
          i1 = pos;
        end
        pos = floor((i1+i2)/2);
      end
    end

  end
end
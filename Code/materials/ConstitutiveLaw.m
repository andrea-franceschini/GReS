classdef ConstitutiveLaw < handle
  %CONSTITUTIVELAW Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    status = struct('curr',[],'conv',[])      % internal state variables
    loc2gp                                    % map local cell index to gp location 
    loc2glob                                  % map region local cell index to global index
  end

  properties (Access=protected)
    ngp
    ncells
  end
  
  methods (Abstract)

    initializeStatus(obj,domain,gaussId)

    constitutiveUpdate(obj)

  end


  methods

    function initialize(obj,tag,domain)

      g = domain.grid;

      % set the cell to gp map
      cells = g.cells;

      cellsInRegion = ismember(cells.tag,tag);

      obj.loc2glob = find(cellsInRegion);

      obj.loc2gp = zeros(length(obj.loc2glob),1);
      obj.loc2gp(:,2) = domain.gpMap(obj.loc2glob,2);
      obj.loc2gp(:,1) = [1;1+cumsum(obj.loc2gp(1:end-1,2))];

      obj.ncells = length(obj.loc2glob);

      % call abstract method of subclasses
      gaussId = getGaussPointIds(domain.gpMap,obj.loc2glob);

      obj.ngp = length(gaussId);

      obj.initializeStatus(domain,gaussId);

      obj.status.conv = obj.status.curr;

    end


    function advanceStatus(obj)

      obj.status.conv = obj.status.curr;

    end


    function goBackStatus(obj)

      obj.status.curr = obj.status.conv;

    end

  end


  methods (Access=protected)

    function gps = getGaussPointIdsFromCell(obj,cellId)

      gps = getGaussPointIds(obj.loc2gp,cellId);
    end

  end
end


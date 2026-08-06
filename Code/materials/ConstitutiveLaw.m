classdef ConstitutiveLaw < handle
  
  properties
    status = struct('curr',[],'conv',[])      % internal state variables
    loc2gp = [1,1]                            % map local cell index to gp location 
  end

  properties (Access=protected)
    ngp = 1
    ncells = 1
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

      loc2glob = find(cellsInRegion);

      obj.loc2gp = zeros(length(loc2glob),1);
      obj.loc2gp(:,2) = domain.gpMap(loc2glob,2);
      obj.loc2gp(:,1) = [1;1+cumsum(obj.loc2gp(1:end-1,2))];

      obj.ncells = length(loc2glob);

      % call abstract method of subclasses
      gaussId = getGaussPointIds(domain.gpMap,loc2glob);

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


  methods (Static)


    function out = isLinear()
      out = false;
    end

  end
end


classdef PeacemanWells < handle
  % Peaceman well models

  % Currently implemented for TPFA flow discretization

  properties (Access = protected)

    wells     % struct of wells defined using MRST addWells
    cellList (:,2)
    domain 

  end

  methods (Access = public)

    
    function obj = PeacemanWells(domain)

      obj.domain = domain;

    end


    function addWell(obj,cells,varargin)

      wId = numel(obj.wells) + 1;

      obj.cellList = [obj.cellList; [repmat(wId,numel(cells),1),cells]];


      % define default structure
      
      newWell = readInput(default,varargin);

      newWell.effRadius = computeEffectiveRadius(obj,newWell);

      newWell.wellIndex = computeWellIndex(obj,newWell);

      obj.wells = [obj.wells; newWell];

    end


    function assembleWells(obj)

      % assemble well contribution within linear system

      
    end


    function initialize(obj)

      % setup the well ADI equations



    end




  end


  methods (Access = private)



    function computeEffectiveRadius(obj,w)
      % ref: White et al. 2019

      dirs = find(strcmp(["x","y","z"]),w.dir);

      if size(dirs) ~= 2
        error("Incorrect definition of well direction")
      end





    


    end


  end



end

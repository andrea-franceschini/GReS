classdef GrowingBoundary < handle
  % STRUCTUREDGRID Class to control all information related to the
  % structured grid

  properties (Access = public)
    growtype    % type of how the boundary it will grow (ONLY Constant)

    entities    % Indices of constrained degrees of freedom in the update BC
    refloc      % Indices of constrained degrees of freedom in relation to the initial mesh
    weight      % Weight of each constrained degrees of freedom in relation to the initial mesh
    numEntBefGw % number of entities before the mesh starts to grow
  end

  properties (Access = private)
    
  end

  methods(Access = public)
    function obj = GrowingBoundary(data,growtype)
      %ADAPGRID Construct an instance of this class
      obj.entities = data.entities;
      obj.numEntBefGw = data.nEntities;
      switch growtype
        case {'null','Null','NULL'}
          obj.growtype = 0;
          obj.weight =[];
        case {'cons','constant','Constant','CONSTANT'}
          obj.growtype = 1;
          obj.refloc = zeros(obj.numEntBefGw,1);
          obj.refloc(:,1)=1:obj.numEntBefGw;
          obj.weight = ones(obj.numEntBefGw,1);
      end
    end

    function update(obj,newCell,refCell,dir)
      %UPDATEBOUNDARY update the boundary condition in respect of the cell

      % Some important shortcut
      ref.dir  = dir;
      ref.dirM = ref.dir+2*mod(ref.dir,2)-1;
      ref.face = newCell.faces(ref.dir);
      ref.faceM = newCell.faces(ref.dirM);
      ref.cellRef = newCell.cell(ref.dir);
      ref.cellNew = refCell.cell(ref.dirM);

      % Start by find the neighbors
      nullCellNeig.New = newCell.cell == 0;
      nullCellNeig.New(ref.dirM) = false;
      nullCellNeig.Ref = refCell.cell == 0;

      isFaceIn = ismember(obj.entities,ref.face);
      if sum(isFaceIn)
        % The reference face is present in this boundary condition.
        loc = find(isFaceIn);

        % Add DOFs in the border
        refLoc = nullCellNeig.New & ~nullCellNeig.Ref;
        nfaces2add = sum(refLoc);
        addFaces = newCell.faces(refLoc)';
        addWeight(1:nfaces2add) = 1;
        addLoc(1:nfaces2add) = obj.refloc(loc);

        obj.entities(loc)=ref.faceM;

        obj.entities(end+1:end+nfaces2add)=addFaces;
        obj.refloc(end+1:end+nfaces2add)=addLoc;
        obj.weight(end+1:end+nfaces2add)=addWeight;

        % Delete DOFs in the border
        for i=1:6
          if i==ref.dir || i==ref.dirM
            continue
          end
          if ~nullCellNeig.New(i)
            loc = find(ismember(obj.entities,newCell.faces(i)));
            obj.entities(loc)=[];
            obj.refloc(loc)=[];
            obj.weight(loc)=[];
          end
        end
      else
        % The reference face is not present in this boundary condition.
        refLoc = nullCellNeig.New & nullCellNeig.Ref;
        refFaces = refCell.faces(refLoc)';
        newFaces = newCell.faces(refLoc)';      
        for i=1:length(refFaces)
          isInside = ismember(obj.entities,refFaces(i));
          if ~sum(isInside)
            continue
          end
          obj.entities(end+1)=newFaces(i);
          obj.refloc(end+1)=obj.refloc(isInside);
          obj.weight(end+1)=obj.weight(isInside);
        end
      end
    end

    function updateBorder(obj,bc)
      nrep = size(obj.refloc,2);
      nent = length(obj.entities);
      availVals=zeros(nent,bc.data.nTimes+1);
      for i=1:2
        vals(1:nent,1) = 0.;
        valsBc = bc.data.availVals(:,i);
        for j=1:nrep
          vals = vals+obj.weight.*valsBc(obj.refloc(:,j));
        end
        availVals(:,i)=vals;
      end
      bc.data.availVals=availVals;
      bc.data.entities=obj.entities;
      bc.data.totEnts=nent;
      bc.data.nEntities=nent;
    end

  end


  methods (Access = private)

  end

  methods (Static)
    
  end

end
classdef Block < handle
  
  properties
    localIndex % local index w.r.t parent 
    globalIndex % global index w.r.t root block
    level
    bbox   % [xmin xmax; ymin ymax; zmin zmax]
    children = Block.empty
    parent
    faceTag = -1*ones(6,1)

  end

  properties (Constant)
    faceToNode = [1,2,3,4
                  1,4,8,5
                  1,2,6,5
                  2,3,7,6
                  4,3,7,8
                  5,6,7,8];
  end

        %        v
    % 4----------3
    % |\     ^   |\
    % | \    |   | \
    % |  \   |   |  \
    % |   8------+---7
    % |   |  +-- |-- | -> u
    % 1---+---\--2   |
    %  \  |    \  \  |
    %   \ |     \  \ |
    %    \|      w  \|
    %     5----------6

  methods
    function obj = Block(bbox,level,localIndex)
      if nargin > 0
        obj.bbox = bbox;
        obj.level = level;
        obj.localIndex = localIndex;
        obj.globalIndex = localIndex;
      end
    end

    function refine(obj)
      % Split current block and build 8 children
      if isempty(obj.children)
        [x, y, z] = obj.subdivideBBox();
        k = 0;
        for l = 1:2
          for j = 1:2
            for i = 1:2
              k = k + 1;
              bb = [x(i:i+1); y(j:j+1); z(l:l+1)];
              locIndex = [i-1 j-1 l-1];
              obj.children(k) = Block(bb, obj.level + 1,locIndex);
              obj.children(k).globalIndex = 2*obj.globalIndex + locIndex;
              obj.children(k).parent = obj;
            end
          end
        end
      end
    end


    function refineBlock(obj,maxDepth,varargin)
      if nargin > 2
        targetPoint = varargin{1};
        if obj.contains(targetPoint)
          if obj.level < maxDepth
            obj.refine();
            for child = obj.children
              refineBlock(child,maxDepth,targetPoint);
            end
          end
        end
      elseif nargin == 2
        if obj.level < maxDepth
          obj.refine();
          for child = obj.children
            refineBlock(child,maxDepth);
          end
        end
      end
    end

    function [x, y, z] = subdivideBBox(obj)
      b = obj.bbox;
      x = linspace(b(1,1), b(1,2), 3);
      y = linspace(b(2,1), b(2,2), 3);
      z = linspace(b(3,1), b(3,2), 3);
    end

    function inside = contains(obj, point)
      b = obj.bbox;
      inside = all([point' >= b(:,1)  point' <= b(:,2)],"all");
    end


    function leaves = getLeaves(obj)
      if isempty(obj.children)
        leaves = obj;
      else
        leaves = [];
        for i = 1:length(obj.children)
          leaves = [leaves, obj.children(i).getLeaves()];  % concatenate
        end
      end
    end

    function coords = getNodeCoordinates(obj)
      % return 8x3 coordinates of nodes based on bounding box coordinates

      coords = [obj.bbox(1,1) obj.bbox(2,1) obj.bbox(3,1);
                obj.bbox(1,2) obj.bbox(2,1) obj.bbox(3,1);
                obj.bbox(1,2) obj.bbox(2,2) obj.bbox(3,1);
                obj.bbox(1,1) obj.bbox(2,2) obj.bbox(3,1);
                obj.bbox(1,1) obj.bbox(2,1) obj.bbox(3,2);
                obj.bbox(1,2) obj.bbox(2,1) obj.bbox(3,2);
                obj.bbox(1,2) obj.bbox(2,2) obj.bbox(3,2);
                obj.bbox(1,1) obj.bbox(2,2) obj.bbox(3,2)];
      
    end


  end
end
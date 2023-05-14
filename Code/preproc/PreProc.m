classdef PreProc < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    indB
    nE   % nE = [#tetra, #hexa, #wed, #pyr]
    nNodesElem = [4, 8, 6, 5]
    nEntryKLoc
    nMat
  end
  
  properties (Access = private)
    material
    GaussPts
    mesh
  end
  
  methods (Access = public)
    function obj = PreProc(grid,mat,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      nIn = nargin;
      data = varargin;
      obj.setPreProc(nIn,grid,mat,data);
    end
    
    function dof = getDoFID(obj,el)
      top = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
      dof = repelem(obj.mesh.nDim*top,obj.mesh.nDim);
      dof = dof + repmat(-(obj.mesh.nDim-1):0,1,obj.mesh.cellNumVerts(el));
    end
    
    function D = getStiffMatrix(obj,el,sz)
      mat = obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw;
      % Material stiffness matrix
      if isa(mat,'HypoElastic')
        D = mat.getStiffnessMatrix(sz);
      else
        D = mat.getStiffnessMatrix();
      end
    end
    
%     function specGrav = getFluidSpecGrav(obj)
%       mat = obj.material.getMaterial(max(obj.mesh.cellTag(el))+1);
%       specGrav = mat.getSpecGrav();
%     end
    
%     function mu = getDynViscosity(obj)
%       mat = obj.material.getMaterial(max(obj.mesh.cellTag(el))+1);
%       mu = mat.getViscosity();
%     end
  end
  
  methods (Access = private)
    function setPreProc(obj,nIn,grid,mat,data)
      obj.mesh = grid.topology;
      obj.material = mat;
      if nIn > 2
        obj.GaussPts = data{1};
      end
      obj.nE = histc(obj.mesh.cellVTKType,[10, 12, 13, 14]);
      obj.nEntryKLoc = (obj.mesh.nDim^2)*(obj.nNodesElem).^2;
      if obj.nE(2) == 0
        l1 = 4;
      else
        l1 = 8*obj.GaussPts.nNode;
      end
      obj.indB = zeros(9*l1,2);
      obj.indB(:,1) = repmat([1, 2, 3, 2, 1, 3, 3, 2, 1],[1,l1]);
      obj.indB(:,2) = repmat([1, 4, 6, 8,10,11,15,17,18],[1,l1]);
      obj.indB(:,1) = obj.indB(:,1) + repelem(3*(0:(l1-1))',9);
      obj.indB(:,2) = obj.indB(:,2) + repelem(18*(0:(l1-1))',9);
      %
      obj.nMat = max(obj.mesh.cellTag);
    end
  end
end
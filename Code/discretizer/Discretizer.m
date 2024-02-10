classdef Discretizer < handle           
  % General discretizer class
  % Database for all models activated in the simulation
  properties (Access = private)
      nField
  end

  properties (Access = public)
    J
    rhs
    fields
    db
  end
  
  properties (Access = public)
    model
    simParams
    dofm
    mesh
    elements
    faces
    material
    GaussPts
    upElem
    fieldMap
  end
  
  methods (Access = public)
      function obj = Discretizer(symmod,simParams,dofManager,grid,mat,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.db = containers.Map('KeyType','char','ValueType','any');
      obj.setDiscretizer(symmod,simParams,dofManager,grid,mat,varargin);
    end
    
        
    function computeVSFMatricesAndRHS(obj,statek,stateTmp,dt)
      pkpt = obj.simParams.theta*stateTmp.pressure + ...
        (1 - obj.simParams.theta)*statek.pressure;
      [Swkpt,dSwkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt);
      computeFlowStiffMatFV(obj,lwkpt);
      computeFlowCapMatFV(obj,Swkpt,dSwkpt);
      computeFlowJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,dlwkpt);
      computeFlowRHS(obj,statek,stateTmp,dt,lwkpt);
    end
    
    function computeFlowJacobian(obj,dt,varargin)
      % varargin -> statek,stateTmp,pkpt,dSwkpt,dlwkpt
      % IF SINGLE PHASE
      obj.J = obj.simParams.theta*obj.H + obj.P/dt;
      if obj.model.isVariabSatFlow() && obj.simParams.isNewtonNLSolver()
        % Compute the additional terms of the Jacobian
        [JNewt] = computeNewtPartOfJacobian(obj,dt,varargin{1}, ...
          varargin{2},varargin{3},varargin{4},varargin{5});
        obj.J = obj.J + JNewt;
      end
    end

    function computeFlowJacobian_Test(obj,dt,varargin)
      nBlock = size(obj.blockJ,1);
          for i=1:nBlock^2
              if all(strcmp(obj.blockJ(i).physics,'Flow')) 
                  % varargin -> statek,stateTmp,pkpt,dSwkpt,dlwkpt
                  % IF SINGLE PHASE
                  obj.blockJ(i).block = obj.simParams.theta*obj.blockJ(i).H + obj.blockJ(i).P/dt;
                  if obj.model.isVariabSatFlow() && obj.simParams.isNewtonNLSolver()
                    % Compute the additional terms of the Jacobian
                    [JNewt] = computeNewtPartOfJacobian(obj,dt,varargin{1}, ...
                      varargin{2},varargin{3},varargin{4},varargin{5});
                    obj.J = obj.J + JNewt;
                  end
              end
          end
    end


    function computeBlockJacobianAndRhs(obj, dt)
        % compute blocks of Jacobian matrix
        for i = 1:obj.nField
            for j = 1:obj.nField
                obj.J{i,j} = obj.getField(obj.fieldMap{i,j}).blockJacobian(i,j,dt);
            end
        end
        % compute blocks of Rhs
        flds = unique(obj.fieldMap);
        for i = 1:obj.nField
            for j = 1:length(flds)
                obj.rhs{i} = obj.rhs{i} + obj.getField(flds{j}).blockRhs(i);
            end
        end

    end

    function resetJacobianAndRhs(obj)
        %reset Jacobian and Rhs to zero
        obj.J = cell(obj.nField,obj.nField);
        obj.rhs = cell(obj.nField,1);
        % set rhs block to zero arrays
        for i = 1:obj.nField
            obj.rhs{i} = zeros(obj.dofm.numDof(i),1);
        end
    end  

    function dSol = solve(obj)
        dSol = cell2mat(obj.J)\(-cell2mat(obj.rhs));
    end

    function out = getSPFlow(obj)
        out = obj.db('Flow');
    end

    function out = getPoro(obj)
        out = obj.db('Poro');
    end

    function out = getBiot(obj)
        out = obj.db('Biot');
    end

  function out = getField(obj,fld)
      out = obj.db(fld);
  end
  end
  
  methods(Access = private)
      function setDiscretizer(obj,symmod,params,dofManager,grid,mat,data)
          obj.model = symmod;
          obj.simParams = params;
          obj.dofm = dofManager;
          obj.material = mat;
          modMap = getModelMap(obj);
          obj.nField = length(dofManager.numDof);
          fldMap = cell(obj.nField, obj.nField);
          for i = 1:obj.nField
              for j = 1:obj.nField
                  str = convertStringsToChars(strcat(dofManager.subPhysics(i),dofManager.subPhysics(j)));
                  fldMap{i,j} = modMap(str);
              end
          end
          obj.fieldMap = fldMap;
          obj.fields = unique(fldMap);
          % initialize block Jacobian and rhs
          obj.J = cell(obj.nField, obj.nField);
          obj.rhs = cell(obj.nField,1);
          % set rhs block to zero arrays
          for i = 1:obj.nField
              obj.rhs{i} = zeros(obj.dofm.numDof(i),1);
          end
          % build model subclasses
          mods = convertCharsToStrings(unique(obj.fieldMap));
          for i = 1:length(mods)
              switch mods(i)
                  case 'Poro'
                      obj.db('Poro') = Poromechanics(symmod,params,dofManager,grid,mat,data);
                  case 'Flow'
                      obj.db('Flow') = SPFlow(symmod,params,dofManager,grid,mat,data);
                  case 'Biot'
                      obj.db('Biot') = Biot(symmod,params,dofManager,grid,mat,data);
              end
          end
          
    
      function mp = getModelMap(obj)
          mp = containers.Map('KeyType','char','ValueType','any');
          mp('PoroPoro') = 'Poro';
          mp('PoroFlow') = 'Biot';
          mp('FlowPoro') = 'Biot';
          mp('FlowFlow') = 'Flow';
      end

    function [JNewt] = computeNewtPartOfJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,dlwkpt)
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      zVec = obj.elements.cellCentroid(:,3);
      zNeigh = zVec(neigh);
      gamma = obj.material.getFluid().getFluidSpecWeight();
      tmpVec1 = (dlwkpt.*obj.trans(obj.isIntFaces)).*(pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)));
      %
      poroMat = zeros(obj.mesh.nCellTag,1);
      alphaMat = zeros(obj.mesh.nCellTag,1);
      beta = obj.material.getFluid().getFluidCompressibility();
      for m = 1:obj.mesh.nCellTag
        poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
        alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
      end
      % (alpha+poro*beta)
      tmpVec2 = alphaMat(obj.mesh.cellTag) + beta*poroMat(obj.mesh.cellTag);
      tmpVec2 = 1/dt*((tmpVec2.*dSwkpt).*(stateTmp.pressure - statek.pressure)).*obj.elements.vol;
%       tmpVec = lw.*tmpVec;
%       sumDiagTrans = accumarray([neigh1; neigh2], ...
%         repmat(tmpVec,[2,1]),[obj.mesh.nCells,1]);
      JNewt = sparse([neigh(:,1); neigh(:,2); (1:obj.mesh.nCells)'], ...
                     [repmat(obj.upElem,[2,1]); (1:obj.mesh.nCells)'], ...
                     [tmpVec1; -tmpVec1; tmpVec2],obj.mesh.nCells,obj.mesh.nCells);
      JNewt = obj.simParams.theta*JNewt;
    end
 
             
    function [Swkpt,dSwkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt)
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      gamma = obj.material.getFluid().getFluidSpecWeight();
      if gamma > 0
        zVec = obj.elements.cellCentroid(:,3);
        zNeigh = zVec(neigh);
        lElemIsUp = pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)) >= 0;
      else
        lElemIsUp = pkpt(neigh(:,1)) >= pkpt(neigh(:,2));
      end
      obj.upElem(lElemIsUp) = neigh(lElemIsUp,1);
      obj.upElem(~lElemIsUp) = neigh(~lElemIsUp,2);
      [Swkpt,dSwkpt] = obj.material.computeSwAnddSw(obj.mesh,pkpt);
%       dSwkpt = zeros(length(dSwkpt),1);
      dSwkpt = - dSwkpt;
      [lwkpt,dlwkpt] = obj.material.computeLwAnddLw(obj.mesh,obj.upElem,pkpt);
      dlwkpt = - dlwkpt;
    end
  end

  end
end
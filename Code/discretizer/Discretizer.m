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
        material
        fieldMap
    end

    methods (Access = public)
        function obj = Discretizer(symmod,simParams,dofManager,grid,mat,varargin)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.db = containers.Map('KeyType','char','ValueType','any');
            obj.setDiscretizer(symmod,simParams,dofManager,grid,mat,varargin);
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

        function out = getFlow(obj)
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

        function computeLinearMatrices(obj,stateTmp,statek,dt)
            for i = 1:length(obj.fields)
                fld = obj.fields{i};
                if isLinear(obj.db(fld))
                    obj.getField(fld).computeMat(stateTmp,statek,dt);
                end
            end
        end

        function computeNLMatricesAndRhs(obj,stateTmp,statek,dt)
            for i = 1:length(obj.fields)
                fld = obj.fields{i};
                if ~isLinear(obj.db(fld))
                    obj.getField(fld).computeMat(stateTmp,statek,dt);
                end
                obj.getField(fld).computeRhs(stateTmp,statek,dt);
            end
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
                  ph1 = translatePhysic(dofManager.subPhysics(i));
                  ph2 = translatePhysic(dofManager.subPhysics(j));
                  str = convertStringsToChars(strcat(ph1,ph2));
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
                      if isSinglePhaseFlow(obj.model)
                          obj.db('Flow') = SPFlow(symmod,params,dofManager,grid,mat,data);
                      end
                      if isVariabSatFlow(obj.model)
                          obj.db('Flow') = VSFlow(symmod,params,dofManager,grid,mat,data);
                      end
                  case 'Biot'
                      if isSinglePhaseFlow(obj.model)
                          obj.db('Biot') = Biot(symmod,params,dofManager,grid,mat,data);
                      end
                      if isVariabSatFlow(obj.model)
                          obj.db('Biot') = BiotVS(symmod,params,dofManager,grid,mat,data);
                      end
              end
          end
          
    
      function mp = getModelMap(obj)
          mp = containers.Map('KeyType','char','ValueType','any');
          mp('PoroPoro') = 'Poro';
          mp('PoroFlow') = 'Biot';
          mp('FlowPoro') = 'Biot';
          mp('FlowFlow') = 'Flow';
      end
  end

  end
end
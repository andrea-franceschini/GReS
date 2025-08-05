classdef FVM < SinglePhysics
   %FVM
   % Subclass of Discretizer for a generic class of Finite Volume Method
   % Implement a generic Finite Volume Methods to assemble the jacobian
   % matrix and the residual contribution for a equation compose by a
   % gradient and/or a divergent.
   %
   % System of the type:
   % $$ Hp + P \frac{\partial p}{\partial t} + f = 0 $$

   properties
      H
      P
      OpDiv     % Matrix
      OpGrad
   end

   methods (Access = public)
      function obj = FVM(symmod,params,dofManager,grid,mat,data,varargin)
         %FVM Construct an instance of this class
         if isempty(varargin)
            fieldName = 'SPFlow';
         else
            fieldName = varargin{1};
         end
         obj@SinglePhysics(fieldName,symmod,params,dofManager,grid,mat,data);
         obj.computeTrans;
         %get cells with active flow model
         flowCells = obj.dofm.getFieldCells(obj.field);
         % Find internal faces (i.e. shared by two cells with active
         % flow)
         obj.isIntFaces = all(ismember(obj.faces.faceNeighbors, flowCells), 2);
         computeRHSGravTerm(obj);
      end

      function [dof,vals] = getBC(obj,bc,id,t,state)
         %GETBC - function to find the value and the location for the
         % boundary condition.
         switch bc.getCond(id)
            case {'NodeBC','ElementBC'}
               ents = bc.getEntities(id);
               vals = bound.getVals(id,t);
            case 'SurfBC'
               [faceID, faceOrder] = sort(bc.getEntities(id));
               ents = sum(obj.faces.faceNeighbors(faceID,:),2);
               [ents,~,ind] = unique(ents);
               switch bc.getType(id)
                  case 'Neu'
                     v(faceOrder,1) = bc.getVals(id,t);
                     area = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
                     vals = accumarray(ind, area);
                  case 'Dir'
                     v(faceOrder,1) = bc.getVals(id,t);
                     [dirJ,q] = findBC(v,state);
                     if isNewtonNLSolver(obj.simParams)
                        [dirJh,qh] = findBCNewtonNL(v,state);
                        dirJ=dirJ+dirJh;
                        q=q+qh;
                     end
                     vals = [dirJ,accumarray(ind,q)];
               end
            case 'VolumeForce'
               v = bc.getVals(id,t);
               ents = bc.getEntities(id);
               vals = v.*obj.elements.vol(ents);
         end
         dof = obj.dofm.getLocalDoF(ents,obj.field);
      end

      function state = computeMat(obj,state,~,dt)
         %GETBC - function to recompute elementary matrices only if the
         % model is non-linear
         if ~isLinear(obj) || isempty(obj.J)
            mu = obj.material.getFluid().getDynViscosity();
            computeStiffMatFV(obj,1/mu);
            computeCapMatFV(obj);
         end
         if obj.simParams.isTimeDependent
            obj.J = obj.simParams.theta*obj.H + obj.P/dt;
         else
            obj.J = obj.H;
         end
      end
   end
   methods (Access = private)
   end
end
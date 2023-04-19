classdef State < matlab.mixin.Copyable
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    t = 0
    iniStress
    stress
    displ
%     iniAvStress
    pressure
  end
  
  properties (Access = private)
    model
    grid
    mesh
    elements
    material
    GaussPts
    preP
%     avStress
%     avStrain
%     fluidPot
%     density = 2200/10^6   % kg/m3/10^6
  end
  
%   properties (Constant)
%     gravAcc = 9.80665   % m/s2
%   end
  
  methods (Access = public)
    function obj = State(symmod,grid,mat,pre,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      nIn = nargin;
      data = varargin;
      obj.setState(nIn,symmod,grid,mat,pre,data);
      obj.iniState();
    end
    
    function updateState(obj,dSol)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      if isPoromechanics(obj.model)
        obj.displ = obj.displ + dSol;
        %
        % Update stress
        l = 0;
        for el=1:obj.mesh.nCells
          dof = getDoFID(obj.preP,el);
          switch obj.mesh.cellVTKType(el)
            case 10 % Tetra
              N = getDerBasisF(obj.elements,el);
              B = zeros(6,4*obj.mesh.nDim);
              B(obj.preP.indB(1:36,2)) = N(obj.preP.indB(1:36,1));
              D = obj.preP.getStiffMatrix(el,obj.stress(l+1,3) ...
                  + obj.iniStress(l+1,3));  % obj.stress before being updated
              dStress = D*B*dSol(dof);
              obj.stress(l+1,:) = obj.stress(l+1,:) + dStress';
  %             obj.avStress(el,:) = obj.stress(l+1,:);
              l = l + 1;
            case 12 % Hexa
              N = getDerBasisFAndDet(obj.elements.hexa,el,2);
              B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
              B(obj.preP.indB(:,2)) = N(obj.preP.indB(:,1));
              D = obj.preP.getStiffMatrix(el,obj.stress(l+1:l+obj.GaussPts.nNode,3) ...
                  + obj.iniStress(l+1:l+obj.GaussPts.nNode,3));  % obj.stress before being updated
              dStress = pagemtimes(D,pagemtimes(B,dSol(dof)));
              obj.stress((l+1):(l+obj.GaussPts.nNode),:) = ...
                obj.stress((l+1):(l+obj.GaussPts.nNode),:) + ...
                reshape(dStress,6,obj.GaussPts.nNode)';
  %             vol = getVolume(obj.elements,el); % Volumetric average
  %             dStress = dStress.*reshape(dJWeighed,1,1,[]);
  %             obj.avStress(el,:) = sum(dStress,3)/vol;
              l = l + obj.GaussPts.nNode;
          end
        end
      end
      %
      if isSinglePhaseFlow(obj.model)
        obj.pressure = obj.pressure + dSol;
      end
    end
    
    function [avStress,avStrain] = finalizeStatePoro(obj)
      avStress = zeros(obj.mesh.nCells,6);
      avStrain = zeros(obj.mesh.nCells,6);
      l = 0;
      for el=1:obj.mesh.nCells
        dof = getDoFID(obj.preP,el);
        switch obj.mesh.cellVTKType(el)
          case 10 % Tetra
            N = getDerBasisF(obj.elements,el);
            B = zeros(6,4*obj.mesh.nDim);
            B(obj.preP.indB(1:36,2)) = N(obj.preP.indB(1:36,1));
            dStrain = B*obj.displ(dof);
            avStrain(el,:) = dStrain;
            avStress(el,:) = obj.stress(l+1,:);
            l = l + 1;
          case 12 % Hexa
            [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
            B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
            B(obj.preP.indB(:,2)) = N(obj.preP.indB(:,1));
%             vol = getVolume(obj.elements,el); % Volumetric average
            avStress(el,:) = sum(diag(dJWeighed)*obj.stress((l+1):(l+obj.GaussPts.nNode),:))/obj.elements.vol(el);
            dStrain = pagemtimes(B,obj.displ(dof));
            dStrain = dStrain.*reshape(dJWeighed,1,1,[]);
            avStrain(el,:) = sum(dStrain,3)/obj.elements.vol(el);
            l = l + obj.GaussPts.nNode;
        end
      end
    end
    
    function [fluidPot] = finalizeStateFlow(obj)
      fluidPot = obj.pressure;
      gamma = obj.material.getMaterial(2*obj.preP.nMat+1).getFluidSpecWeight();
      if gamma > 0
        if isFEMBased(obj.model)
          fluidPot = fluidPot + gamma*obj.mesh.coordinates(:,3);
        elseif isFVTPFABased(obj.model)
          fluidPot = fluidPot + gamma*obj.elements.cellCentroid(:,3);
        end
      end
    end
    
    function transferState(objFrom,objTo)
      %
      % Shared
      objTo.t = objFrom.t;
      % Poromechanics
      objTo.displ = objFrom.displ;
      objTo.stress = objFrom.stress;
      %
      % Flow
      objTo.pressure = objFrom.pressure;
    end
    
%     function outStress(obj)
%       if sum(obj.preP.nE(2:end)) == 0
%         obj.avStress = obj.stress;
%       else
%         l = 0;
%         for el=1:obj.mesh.nCells
%           switch obj.mesh.cellVTKType(el)
%             case 10 % Tetra
%               obj.avStress(el,:) = obj.stress(l+1,:);
%               l = l + 1;
%             case 12 % Hexa
%               
%               l = l + obj.GaussPts.nNode;
%           end
%         end
%       end
%     end
  end
  
  methods (Access = private)
    function setState(obj,nIn,symmod,grid,mat,pre,data)
      obj.model = symmod;
      obj.grid = grid;
      obj.mesh = grid.topology;
      obj.elements = grid.cells;
      obj.material = mat;
      obj.preP = pre;
%       obj.density = dens;
      if nIn > 4
        obj.GaussPts = data{1};
      end
      if isPoromechanics(obj.model)
        tmp = 1;
        if nIn > 4
          tmp = obj.GaussPts.nNode;
        end
        % NOT ELEGANT AT ALL. FIX!
        obj.stress = zeros([1, tmp, 0, 0]*obj.preP.nE,6);
        obj.iniStress = zeros([1, tmp, 0, 0]*obj.preP.nE,6);
        obj.displ = zeros(obj.mesh.nDim*obj.mesh.nNodes,1);
%         obj.avStress = zeros(obj.mesh.nCells,6);
%         obj.avStrain = zeros(obj.mesh.nCells,6);
%         obj.iniAvStress = zeros(obj.mesh.nCells,6);
      end
      %
      if isSinglePhaseFlow(obj.model)
        if isFEMBased(obj.model,'Flow')
          nr = obj.mesh.nNodes;
        elseif isFVTPFABased(obj.model,'Flow')
          nr = obj.mesh.nCells;
        end
        obj.pressure = zeros(nr,1);
      end
    end
    
    function iniState(obj)
      max_z = max(obj.mesh.coordinates(:,3));
      if isPoromechanics(obj.model)
        l1 = 0;
        for el = 1:obj.mesh.nCells
          M = obj.material.getMaterial(obj.mesh.cellTag(el)).getMFactor();
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          specGrav = 0.0216;    % FIX THE CALL TO THE PROPERTY IN MATERIAL - POROUS ROCK
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          switch obj.mesh.cellVTKType(el)
            case 10 % Tetra
              obj.iniStress(l1+1,3) = -specGrav*(max_z-obj.elements.cellCentroid(el,3));
              obj.iniStress(l1+1,2) = M*obj.iniStress(l1+1,3);
              obj.iniStress(l1+1,1) = obj.iniStress(l1+1,2);
%               obj.iniAvStress(el,:) = obj.iniStress(l1+1,:);
              l1 = l1 + 1;
            case 12 % Hexa
              %%%%%%%%%%%%%%%%%%%%%%%%%% FIX CALL TO TWO ELEMENTS %%%%%%%%%
%               [~,dJWeighed] = getDerBasisFAndDet(obj.elements,el);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              GPLoc = obj.elements.hexa.getGPointsLocation(el);
              obj.iniStress(l1+1:l1+obj.GaussPts.nNode,3) = -specGrav*(max_z-GPLoc(:,3))-50;  %MPa
              obj.iniStress(l1+1:l1+obj.GaussPts.nNode,1) = M*obj.iniStress(l1+1:l1+obj.GaussPts.nNode,3);
              obj.iniStress(l1+1:l1+obj.GaussPts.nNode,2) = obj.iniStress(l1+1:l1+obj.GaussPts.nNode,1);
%               vol = getVolume(obj.elements,el);
%               obj.iniAvStress(el,:) = ((obj.iniStress(l1+1:l1+obj.GaussPts.nNode,:))'*dJWeighed')/vol;
              l1 = l1 + obj.GaussPts.nNode;
          end
        end
      end
      %
      if isSinglePhaseFlow(obj.model)
%         if 5<1
        max_z = 410;
        gamma = obj.material.getMaterial(2*obj.preP.nMat+1).getFluidSpecWeight();
        if isFEMBased(obj.model,'Flow')
%           obj.pressure = gamma*(max_z-obj.mesh.coordinates(:,3));
        obj.pressure = zeros(length(obj.mesh.coordinates(:,3)),1);
%         obj.pressure = 392.4*ones(length(obj.mesh.coordinates(:,3)),1);
        elseif isFVTPFABased(obj.model,'Flow')
%           obj.pressure = gamma*(max_z-obj.elements.cellCentroid(:,3));
          obj.pressure = zeros(obj.mesh.nCells,1);
        end
%         end
      end
    end
  end
  
  methods (Access = protected)
    function cp = copyElement(obj)
      if isempty(obj.GaussPts)
        cp = State(obj.model,obj.grid,obj.material,obj.preP);
      else
        cp = State(obj.model,obj.grid,obj.material,obj.preP,obj.GaussPts);
      end
      %
      cp.stress = obj.stress;
      cp.displ = obj.displ;
%       cp.avStress = obj.avStress;
      cp.iniStress = obj.iniStress;
      cp.displ = obj.displ;
%       cp.iniAvStress = obj.iniAvStress;
      cp.pressure = obj.pressure;
      cp.t = obj.t;
    end
  end
end
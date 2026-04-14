classdef (Abstract) MortarQuadrature < handle
  % General class for handling mortar integration
  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    interfacePairs
    grids                  % handle to master and slave lower dimensional grids
    gpCoords = cell(1,2);  % store local coordinates for each master and slave pair
    numbInterfacePairs
    multiplierType
    areaSlave
    gaussOrder
  end

  methods
    function obj = MortarQuadrature(multType,grids,input)
      % grids: the lower dimensional grids for slave and master side

      input = readInput(struct('gaussOrder',3),input);
      obj.gaussOrder = input.gaussOrder;

      obj.multiplierType = multType;
      obj.grids = grids;

    end


  end



  methods (Abstract)
    processMortarPairs(obj,connectivity)
    isPairActive = processMortarPair(obj,is,im);
    finalizeMortarMaps(obj);
    getIntegrationWeights(obj,idPair);
    xiMaster = getMasterGPCoords(obj,idPair);
    xiSlave = getSlaveGPCoords(obj,idPair);
  end

  methods

    % function elem = getElem(obj,sideID,id)
    %   % get instance of element class on one the sides of the interface
    %   % Assumption: same element type in the entire interface
    %   % get istance of element class based on cell type
    %   type = obj.msh(sideID).surfaceVTKType(id);
    %   elem = obj.elements(sideID).getElement(type);
    % end

    function Nmult = computeMultiplierBasisF(obj,elId,elem,NslaveIn)

      switch obj.multiplierType
        case 'P0'
          Nmult = ones(size(NslaveIn,1),1);
        case 'standard'
          Nmult = NslaveIn;
        case 'dual'
          % ref: Popp, A. (2012). Mortar methods for computational contact
          % mechanics and general interface problems 
          Ns = getBasisFinGPoints(elem);
          gpW = getDerBasisFAndDet(elem,elId);
          Ml = Ns'*(Ns.*gpW');
          Dl = diag(Ns'*gpW');
          A = Ml\Dl;
          Nmult = NslaveIn*A;
      end
    end

    function [Nslave, Nmaster, Nmult, varargout] = getMortarBasisFunctions(obj,~,is,elemMaster,elemSlave,xiMaster,xiSlave)

      Nmaster = elemMaster.computeBasisF(xiMaster);
      Nslave = elemSlave.computeBasisF(xiSlave);
      Nmult = obj.computeMultiplierBasisF(is,elemSlave,Nslave);

      if nargout ==4
        % outout the bubble function on the slave side
        varargout{1} = elemSlave.computeBubbleBasisF(xiSlave);
      end
    end

    function computeAreaSlave(obj)
      % compute the effective area of the slave elements being integrated
      obj.areaSlave = zeros(obj.grids(MortarSide.slave).surfaces.num,1);
      for iPair = 1:size(obj.interfacePairs,1)
        is = obj.interfacePairs(iPair,1);
        areaPair = sum(obj.getIntegrationWeights(iPair));
        obj.areaSlave(is) = obj.areaSlave(is) + areaPair;
      end
      
    end


  end

  methods (Static)

    function mat = integrate(func,varargin)
      dJw = varargin{end};
      mat = func(varargin{1:end-1});
      mat = mat.*reshape(dJw,1,1,[]);
      mat = sum(mat,3);
    end
    
  end

    
end


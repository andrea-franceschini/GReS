classdef (Abstract) MortarQuadrature < handle
  % General class for handling mortar integration
  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    msh       % instance of InterfaceMesh class
    interface    % instance of the interface solver object 
    interfacePairs
    gpCoords = cell(1,2);  % store local coordinates for each master and slave pair
    numbInterfacePairs
    elements
    multiplierType
    areaSlave
    
  end

  methods
    function obj = MortarQuadrature(interf,multType,input)

      obj.interface = interf;
      obj.msh = interf.interfMesh.msh;
      obj.multiplierType = multType;

      ng = getXMLData(input,4,"nGP");


      if strcmp(class(obj),"SegmentBasedQuadrature") && ng > 6
        ng = 6;
      end

      obj.elements = [Elements(interf.interfMesh.msh(1),ng),...
        Elements(interf.interfMesh.msh(2),ng)];
    end


  end



  methods (Abstract)
    processMortarPairs(obj)
    isPairActive = processMortarPair(obj,is,im);
    finalizeMortarMaps(obj);
    getIntegrationWeights(obj,idPair);
    xiMaster = getMasterGPCoords(obj,idPair);
    xiSlave = getSlaveGPCoords(obj,idPair);
  end

  methods

    function elem = getElem(obj,sideID,id)
      % get instance of element class on one the sides of the interface
      % Assumption: same element type in the entire interface
      % get istance of element class based on cell type
      type = obj.msh(sideID).surfaceVTKType(id);
      elem = obj.elements(sideID).getElement(type);
    end

    function Nmult = computeMultiplierBasisF(obj,el,NslaveIn)
      elem = obj.getElem(2,el);
      switch obj.multiplierType
        case 'P0'
          Nmult = ones(size(NslaveIn,1),1);
        case 'standard'
          Nmult = NslaveIn;
        case 'dual'
          % ref: Popp, A. (2012). Mortar methods for computational contact
          % mechanics and general interface problems 
          Ns = getBasisFinGPoints(elem);
          gpW = getDerBasisFAndDet(elem,el);
          Ml = Ns'*(Ns.*gpW');
          Dl = diag(Ns'*gpW');
          A = Ml\Dl;
          Nmult = NslaveIn*A;
      end
    end

    function [Nslave, Nmaster, Nmult, varargout] = getMortarBasisFunctions(obj,im,is,xiMaster,xiSlave)
      elemMaster = obj.getElem(1,im);
      elemSlave = obj.getElem(2,is);
      Nmaster = elemMaster.computeBasisF(xiMaster);
      Nslave = elemSlave.computeBasisF(xiSlave);
      Nmult = obj.computeMultiplierBasisF(is,Nslave);

      if nargout ==4
        % outout the bubble function on the slave side
        varargout{1} = elemSlave.computeBubbleBasisF(xiSlave);
      end
    end

    function computeAreaSlave(obj)
      % compute the effective area of the slave elements being integrated
      obj.areaSlave = zeros(obj.msh(2).nSurfaces,1);
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


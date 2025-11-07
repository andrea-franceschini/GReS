classdef (Abstract) MortarQuadrature < handle
  % General class for handling mortar integration
  % REFS: Puso,2004, A mortar segment-to-segment contact method for large
  % deformation solid mechanics
  
  properties
    msh       % instance of InterfaceMesh class
    mortar    % instance of the mortar object 
    mortarPairs
    gpCoords = cell(1,2);  % store local coordinates for each master and slave pair
    numbMortarPairs
    
  end

  methods
    function obj = MortarQuadrature(mortar,ng)
      obj.mortar = mortar;
      obj.msh = mortar.mesh.msh;
      obj.mortar.elements = [Elements(mortar.mesh.msh(1),ng),...
                             Elements(mortar.mesh.msh(2),ng)];
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

  methods (Static)

    function mat = integrate(func,varargin)
      dJw = varargin{end};
      mat = func(varargin{1:end-1});
      mat = mat.*reshape(dJw,1,1,[]);
      mat = sum(mat,3);
    end
    
  end

    
end


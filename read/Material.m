classdef Material < handle
    %MATERIAL General material class
    
    properties
        
        %Material mathematical model TAG, for the definition of the material's
        %properties
        MatModel = 0;
        
        %Total number of materials of the mesh
        nMat = 0;

        %Material ID, must be different for each material
        MatID = [];

        %General properties:
        E = [];             %Elastic modulus 
        ni = [];            %Poisson ratio 
        permx = [];         %Permeability x 
        permy = [];         %Permeability y 
        permz = [];         %Permeability z 
        por = [];           %Porosity
        gammaw = [];        %Material weight 
        betaf = [];         %Fluid compressibility?
        alpha = [];         %Biot coefficient?
        a = [];             %?
        b = [];             %?
        e0 = [];            %Void index
        Cc = [];            %Compression coefficient
        Cr = [];            %Recompression index
        %OCR = [];          %Over Consolidation Ratio
       
    end
    
    properties (Access = private)
        
        %Material Models' TAG
        %1 = 'lin-elas'
        %2 = 'nonlin-elas'
        
    end

    methods
        function setMaterial(obj,fileName)
            %getting material's properties from INPUT file
            [numMat,MatData] = readMaterial(fileName);
            %assigning number of materials and materials' IDs
            obj.nMat = numMat;
            obj.MatID = MatData(:,1);
            %assigning material's maathematical model
            obj.MatModel = MatData(:,2);
            
            %for each material:
            for i = 1:numMat
                %assigning properties 
                        obj.permx = MatData(:,3);
                        obj.permy = MatData(:,4);
                        obj.permz = MatData(:,5);
                        obj.por = MatData(:,6);
                        obj.betaf = MatData (:,7);
                        obj.gammaw = MatData(:,8);
                        obj.E = MatData(:,9);
                        obj.ni = MatData(:,10);
                        obj.alpha = MatData(:,11);
                        
                      %switching properties for different mathematical model  
                      if obj.MatModel(i) == 1 %lin-elas
                        %no properties to add
                        
                        %non lin properties in progress
                      elseif obj.MatModel(i) == 2 %nonlin-elas
                        obj.a(i) = MatData(i,12);
                        obj.b(i) = MatData(i,13);
                        obj.e0(i) = MatData(i,14);
                        obj.Cc(i) = MatData(i,15);
                        obj.Cr(i) = MatData(i,16);
                        
                      end
                end
            
            end
        end

end

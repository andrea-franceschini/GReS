classdef Materials < handle
  % MATERIAL - General material class

  properties (Access = public)
    % Creation of a Dictionary object (faster than map)
    %db = configureDictionary("double","struct"); 
    % configureDictionary is not supported before 2023b
    db 
    matMap
  end

  methods (Access = public)
    % Class constructor method   
    function obj = Materials(model,fListName)
      obj.db = containers.Map('KeyType','double','ValueType','any');
      % Calling the function to read input data from file
      obj.matMap = zeros(100,100);
      obj.readInputFiles(model,fListName)
    end

    % Get the material defined by matIdentifier and check if it is a
    % key of the Map db
    function mat = getMaterial(obj,cellID)
      %
      % The preliminary check whether matID is key of db has been commented
      % since it is highly expensive
%       if (obj.db.isKey(matID))
        [matID,~] = find(obj.matMap == cellID);
        assert(isscalar(matID),['Zero or Multiple materials assigned to elements',...
            ' with cellTags %i'], cellID)
        mat = obj.db(matID);
%       else
%       Displaying error message if matIdentifier is not a key 
%       of the map
%         error('Material %s not present', matID);
%       end
    end

    function fluidMat = getFluid(obj)
        %fluid materials corresponds to null rows in materials map. 
        %Fluid are not assigned to any cellTag
        f = find(sum(obj.matMap,2)==0);
        fluidMat = obj.db(f);
    end
    
    function varargout = computeSwAnddSw(obj,mesh,pkpt)
        % Consider moving this method to a more appropriate location
        % varargout{1} -> Sw
        % varargout{2} -> dSw
        % varargout{3} -> d2Sw
        % if 5<1
        varargout{1} = zeros(mesh.nCells,1);
        if nargout == 2
            varargout{2} = zeros(mesh.nCells,1);
        end
        if nargout == 3
            varargout{2} = zeros(mesh.nCells,1);
            varargout{3} = zeros(mesh.nCells,1);
        end
        for m = 1:mesh.nCellTag
            isElMat = mesh.cellTag == m;
            p = pkpt(isElMat);
            Sws = obj.getMaterial(m).PorousRock.getMaxSaturation();
            Swr = obj.getMaterial(m).PorousRock.getResidualSaturation();            
            % Sws = 1.;
            % Sws = 0.368;
            if nargout == 1
                [varargout{1}(isElMat), ~, ~] = obj.getMaterial(m).Curves.computeSwAnddSw(p,Swr,Sws);
                % varargout{1}(isElMat) = obj.getMaterial(m).CapillaryCurve.interpTable(p);
                % varargout{1}(isElMat) = Swr + (Sws-Swr)*varargout{1}(isElMat);
            elseif nargout == 2
                [varargout{1}(isElMat), varargout{2}(isElMat), ~] = obj.getMaterial(m).Curves.computeSwAnddSw(p,Swr,Sws);
                % [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).CapillaryCurve.interpTable(p);
                % varargout{1}(isElMat) = Swr + (Sws-Swr)*varargout{1}(isElMat);
                % varargout{2}(isElMat) = (Sws-Swr)*varargout{2}(isElMat);
            elseif nargout == 3
                [varargout{1}(isElMat), varargout{2}(isElMat), varargout{3}(isElMat)] = obj.getMaterial(m).Curves.computeSwAnddSw(p,Swr,Sws);
                % [varargout{1}(isElMat), varargout{2}(isElMat), varargout{3}(isElMat)] = obj.getMaterial(m).CapillaryCurve.interpTable(p);
                % varargout{1}(isElMat) = Swr + (Sws-Swr)*varargout{1}(isElMat);
                % varargout{2}(isElMat) = (Sws-Swr)*varargout{2}(isElMat);
                % varargout{3}(isElMat) = (Sws-Swr)*varargout{3}(isElMat);                

                % [St, dSt, ddSt] = obj.getMaterial(m).CapillaryCurve.interpTable(p);
                % St = Swr + (Sws-Swr)*St;
                % dSt = (Sws-Swr)*dSt;
                % ddSt = (Sws-Swr)*ddSt;
                % [Sa, dSa, ddSa] = obj.getMaterial(m).Curves.computeSwAnddSw(p,Swr,Sws);
                % es(1)=sum(Sa-St);
                % es(2)=sum(dSa-dSt);
                % es(3)=sum(ddSa-ddSt);
                % en(1)=norm(Sa-St);
                % en(2)=norm(dSa-dSt);
                % en(3)=norm(ddSa-ddSt);
                % 
                % xx = linspace(-10,-0.75,400);
                % xx = linspace(-1000/980.66,-75/980.66,400);
                % xx = linspace(-1000,-75,400);
                % % xx = linspace(-100,-1e-5,400);
                % % [St, dSt, ddSt] = obj.getMaterial(m).CapillaryCurve.interpTable(xx);
                % % St = Swr + (Sws-Swr)*St;
                % % dSt = (Sws-Swr)*dSt;
                % % ddSt = (Sws-Swr)*ddSt;
                % [Sa, dSa, ddSa] = obj.getMaterial(m).Curves.computeSwAnddSw(xx,Swr,Sws);
                % figure();
                % hold on;
                % plot(xx,Sa,'b', 'LineWidth', 2, 'MarkerSize', 10);
                % % plot(-xx,ddSt,'r', 'LineWidth', 2, 'MarkerSize', 10);
                % xlabel('Pressure');
                % ylabel('Saturation');
                % legend('Analytical','Tabular')
                % xlim([min(xx),max(xx)]);
                % ylim([Swr,Sws]);
                % grid on
                % % axes('Position',[.35 .65 .2 .2])
                % axes('Position',[.25 .25 .3 .3])
                % box on
                % hold on
                % plot(-xx,ddSa,'b', 'LineWidth', 2, 'MarkerSize', 10);
                % plot(-xx,ddSt,'r', 'LineWidth', 2, 'MarkerSize', 10);
                % % xlim([10,12]);
                % xlim([10,15]);
                % hold off
                % set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
                % a = get(gca,'XTickLabel');
                % set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 12)
                % % xlim([-1000,-75]);
                % % ylim([0.102,0.2]);
                % % xticks([-1000,-800,-600,-400,-200]); % Set x-axis ticks
                % % yticks([0.12,0.14,0.16,0.18,0.2]); % Set y-axis ticks
                % % 
                % % hold off;
                % % figure();
                % % hold on;
                % % plot(-xx,Sa-St','b');
                % % plot(-xx,dSa-dSt','r');
                % % plot(-xx,ddSa-ddSt','g');
                % % xlabel('Pressure');
                % % ylabel('Error');
                % % legend('Saturation','First Derivative','Second Derivative');
                % % hold off;
                % % es(1)=sum(Sa-St');
                % % es(2)=sum(dSa-dSt');
                % % es(3)=sum(ddSa-ddSt');
                % % en(1)=norm(Sa-St');
                % % en(2)=norm(dSa-dSt');
                % % en(3)=norm(ddSa-ddSt');
            end
            
            
        end
    end
    
    function varargout = computeLwAnddLw(obj,mesh,upElem,pkpt)
      % Consider moving this method to a more appropriate location
      % varargout{1} -> lw
      % varargout{2} -> dlw
      nIntFaces = length(upElem);
      varargout{1} = zeros(nIntFaces,1);
      if nargout == 2
        varargout{2} = zeros(nIntFaces,1);
      end
      matUpElem = mesh.cellTag(upElem);
      for m = 1:mesh.nCellTag
        isElMat = matUpElem == m;
        p = pkpt(upElem(isElMat));
        if nargout == 1
          [varargout{1}(isElMat), ~] = obj.getMaterial(m).Curves.computeRelativePermeability(p);
          % varargout{1}(isElMat) = obj.getMaterial(m).RelativePermCurve.interpTable(p);
        elseif nargout == 2
          [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).Curves.computeRelativePermeability(p);
          % [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).RelativePermCurve.interpTable(p);          

          % [ka, dka] = obj.getMaterial(m).Curves.computeRelativePermeability(p);
          % [kt, dkt] = obj.getMaterial(m).RelativePermCurve.interpTable(p);          
          % es(1)=sum(ka-kt);
          % es(2)=sum(dka-dkt);
          % en(1)=norm(ka-kt);
          % en(2)=norm(dka-dkt);
          % 
          % xx = linspace(-10,-0.75,400);
          % % xx = linspace(-1000,-0.1,400);
          % % xx = linspace(-100,-1e-5,400);
          % [ka, dka] = obj.getMaterial(m).Curves.computeRelativePermeability(xx);
          % % [kt, dkt] = obj.getMaterial(m).RelativePermCurve.interpTable(xx); 
          % figure();
          % hold on;
          % plot(xx,ka,'b', 'LineWidth', 2, 'MarkerSize', 10);
          % % plot(-xx,dkt,'r', 'LineWidth', 2, 'MarkerSize', 10);
          % xlabel('Pressure');
          % ylabel('Relative Permeability');
          % legend('Analytical','Tabular');
          % xlim([min(xx),max(xx)]);
          % grid on;
          % hold on;
          % axes('Position',[.25 .5 .3 .3])
          % box on
          % hold on
          % plot(-xx,dka,'b', 'LineWidth', 2, 'MarkerSize', 10);
          % plot(-xx,dkt,'r', 'LineWidth', 2, 'MarkerSize', 10);
          % xlim([10,12]);
          % % xlim([-1000,-75]);
          % % ylim([0,3.e-3]);
          % set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
          % a = get(gca,'XTickLabel');
          % set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 12)
          % hold off;

          % figure();
          % hold on;
          % plot(-xx,ka-kt','b');
          % plot(-xx,dka-dkt','r');
          % xlabel('Pressure');
          % ylabel('Error');
          % legend('Relative Permeability','First Derivative');
          % hold off;
          % es(1)=sum(ka-kt');
          % es(2)=sum(dka-dkt');
          % en(1)=norm(ka-kt');
          % en(2)=norm(dka-dkt');
        end
      end
      if (max(p)>-75)
          % pause;
      end

      mu = obj.getFluid().getDynViscosity();
      varargout{1} = varargout{1}/mu;
      varargout{2} = varargout{2}/mu;
    end
    
    function [status] = initializeStatus(obj,cTag,sigma)
      mat = obj.getMaterial(cTag).ConstLaw;
      [status] = mat.initializeStatus(sigma);
    end
    
    function [D, sigma, status] = updateMaterial(obj, cTag, sigma, epsilon, dt, status, el, t)
       mat = obj.getMaterial(cTag).ConstLaw;
       [D, sigma, status] = mat.getStiffnessMatrix(sigma, epsilon, dt, status, el);
    end

    %
    % Destructor
    function delete(obj)
      remove(obj.db,keys(obj.db));
    end
  end

  methods (Access = private)
    % Reading the input file by material blocks
    function readInputFile(obj, model, matFileName, matID, cellTags)
      fID = Materials.openReadOnlyFile(matFileName);
      %
      assert(~feof(fID),'No material properties have been assigned in %s',matFileName);
      propName = readToken(fID, matFileName);
      while true
        % Calling the specific material class based on the
        % material name
        switch propName
          case 'Elastic'
            matProp.ConstLaw = Elastic(fID, matFileName);
          case 'HypoElastic'
            matProp.ConstLaw = HypoElastic(fID, matFileName);
          case 'TransvElastic'
            matProp.ConstLaw = TransvElastic(fID, matFileName);
          case 'SSCM'
            matProp.ConstLaw = SSCM(fID, matFileName);
          case 'PorousRock'
            matProp.PorousRock = PorousRock(fID, model, matFileName);
          case 'CapillaryCurve'
            matProp.CapillaryCurve = TabularCurve(fID, matFileName);
          case 'RelativePermCurve'
            matProp.RelativePermCurve = TabularCurve(fID, matFileName);
          case 'VanGenuchten'
            matProp.Curves = VanGenuchten(fID, matFileName);
          case 'Fluid'
            matProp = Fluid(fID, matFileName);
          otherwise
            error('Error in file %s\nMaterial property %s not available', matFileName, propName);
        end
        token = readToken(fID, matFileName);
        if isempty(strtrim(token))
          token = readToken(fID, matFileName);
          assert(~isempty(token),'Multiple newlines after property %s\nIn material file %s', propName, matFileName);
          propName = token;
        elseif strcmp(token,'End')
          break
        end
      end
      if isfield(matProp,'ConstLaw')
          %check if non-fluid material has cellTags assigned
          assert(~isempty(cellTags),"Undefined cellTags for material #%i", matID)
      else
          assert(isempty(cellTags),"Fluid material do not require any cellTag", matID)
      end
      obj.db(matID) = matProp;
      fclose(fID);
    end
         
    % Reading boundary input file
    function readInputFiles(obj, model, fListName)
      fID = Materials.openReadOnlyFile(fListName);
      [cellTags,matFileName] = readTokenList(fID, fListName);
      matID = 0;
      while ~strcmp(matFileName,'End')
        matID = matID + 1;
        readInputFile(obj, model, matFileName, matID, cellTags);
        obj.matMap(matID,1:length(cellTags)) = cellTags;
        [cellTags,matFileName] = readTokenList(fID, fListName);
      end
      obj.matMap = obj.matMap(1:matID,sum(obj.matMap,1)~=0);
      fclose(fID);

      % Test necessity of beta correction in Van Genuchten object.
      if obj.getFluid().isvalid || obj.getMaterial(1).Curves.isvalid
          if obj.getMaterial(1).Curves.needBetaCor()
            specWeight = obj.getFluid().getFluidSpecWeight();
            obj.db(1).Curves.betaCorrection(specWeight);
          end
      end
    end
  end

  methods (Static = true)
    % Open a read-only file and check  
    function fID = openReadOnlyFile(fName)
      fID = fopen(fName, 'r');
      if fID == -1
        error('File %s does not exist in the directory',fName);
      end
    end

  end
end
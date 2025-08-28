classdef VGMCurves < handle
  %UNTITLED2 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = private)
    soilVGMTable
    soilVGMParams
  end
  
  methods (Access = public)
    function obj = VGMCurves()
      setSoilDatabase(obj);
    end
    
    function setVGMParameters(obj,name,varargin)
      assert(rem(length(varargin), 2) == 0, 'Unformatted input in setVGMParameters');
      if ismember(name,obj.soilVGMTable.Properties.RowNames)
        tmpVec = obj.soilVGMTable{name,:};
        obj.soilVGMParams.n = tmpVec(1);
        obj.soilVGMParams.pEntry = ...
          tmpVec(2)*VGMCurves.setOptions(varargin, 'SpecWeight', 'scalar', 1);
      elseif strcmp(name,'User_defined')
        % Set the parameters according to the input
        obj.soilVGMParams = VGMCurves.procInputArgVGMParam(varargin);
      else
        error('Unknown soil type');
      end
    end
    
    function makeCurves(obj,varargin)
      discrInput = VGMCurves.procInputArgVGMCurve(varargin);
      m = 1 - 1/obj.soilVGMParams.n;
      SeFun = @(p) (1 + (p/obj.soilVGMParams.pEntry).^obj.soilVGMParams.n).^(-m);
      krFun = @(p) (1 + (p/obj.soilVGMParams.pEntry).^obj.soilVGMParams.n).^(-2.5.*m) .* ...
        ((1 + (p/obj.soilVGMParams.pEntry).^obj.soilVGMParams.n).^m - ...
        ((p/obj.soilVGMParams.pEntry).^obj.soilVGMParams.n).^m).^2;
      % points = logspace(discrInput.range(1),discrInput.range(2),discrInput.nPoints);
      points = linspace(10^discrInput.range(1),10^discrInput.range(2),discrInput.nPoints);
      points = [0 points];
      SePoints = SeFun(points);
      krPoints = krFun(points);
      %
      figure('Position', [100, 100, 700, 700])
      fplot(SeFun,[10^(discrInput.range(1)),10^(discrInput.range(2))],'r-');
      %set(gca,'XScale','log');
      hold on
      plot(points(2:end),SePoints(2:end),'r*', 'MarkerSize', 14);
      title('Capillary curve');
      xlabel('Pressure');
      ylabel('S_e');
      xlim([10^(discrInput.range(1)),10^(discrInput.range(2))]);
      ylim([-0.01, 1.01]);
      set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
      hold off
      %
      figure('Position', [100, 100, 700, 700])
      fplot(krFun,[10^(discrInput.range(1)),10^(discrInput.range(2))],'r-');
      %set(gca,'XScale','log','YScale','log');
      hold on
      plot(points(2:end),krPoints(2:end),'r*', 'MarkerSize', 14);
      title('Relative permeability curve');
      xlabel('Pressure');
      ylabel('k_r');
      xlim([10^(discrInput.range(1)),10^(discrInput.range(2))]);
      ylim([0, 1]);
      set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
      hold off
      %
      % Print curves
      fIdPC = fopen(discrInput.fileName(1),'w');
      if fIdPC == -1
        error('File %s not opened correctly',discrInput.fileName(1));
      end
      %
      fIdkr = fopen(discrInput.fileName(2),'w');
      if fIdkr == -1
        error('File %s not opened correctly',discrInput.fileName(2));
      end
      %
      for i=1:length(points)
        fprintf(fIdPC,'%25.15e    %25.15e\n',points(i),SePoints(i));
        fprintf(fIdkr,'%25.15e    %25.15e\n',points(i),krPoints(i));
      end
      fclose(fIdPC);
      fclose(fIdkr);
    end
  end
  
  methods (Access = private)
    function setSoilDatabase(obj)
      names = ["Clay"; "Clayey_loam"; "Loam"; "Loamy_sand"; "Sand"; ...
        "Sandy_clay"; "Sandy_clayey_loam"; "Sandy_loam"; "Silt"; ...
        "Silty_clay"; "Silty_clayey_loam"; "Silty_loam"];
      n = [0.098; 0.151; 0.168; 0.242; 0.502; 0.082; 0.124; 0.161; ...
        0.225; 0.121; 0.182; 0.221];
      n = 10.^(n);
      alpha = [-1.825; -1.801; -1.954; -1.459; -1.453; -1.476; -1.676; ...
        -1.574; -2.182; -1.790; -2.076; -2.296];
      alpha = 10.^(-alpha)/100;
      obj.soilVGMTable = table(n,alpha,'RowNames',names);
    end
  end
  
  methods (Static = true)
    function opt = setOptions(data, name, varType, len)
      id = strcmp((string(data(1:2:end))),name);
      if sum(id) ~= 1
        error('%s not defined or defined multiple times', name);
      end
      switch varType
        case 'string'
          assert(all(isstring(data{2*find(id)})),'%s is not a string',name);
        case 'number'
          assert(all(isnumeric(data{2*find(id)})),'%s is not a number',name);
      end
      assert(length(data{2*find(id)}) == len, 'Unexpected number of entries in %s', name);
      opt = data{2*find(id)};
    end
    
    function out = procInputArgVGMCurve(data)
      assert(rem(length(data),2) == 0,'Unformatted input in makeCurves');
      out.fileName = VGMCurves.setOptions(data, 'fileName', 'string', 2);
      out.range = VGMCurves.setOptions(data, 'range', 'number', 2);
      out.nPoints = VGMCurves.setOptions(data, 'nPoints', 'number', 1);
    end
    
    function out = procInputArgVGMParam(data)
      out.n = VGMCurves.setOptions(data, 'n', 'number', 1);
      if out.n < 1.25 || out.n > 6
        warning('The n value is outside the suggested range (1.25 < n < 6)');
      end
      out.pEntry = VGMCurves.setOptions(data, 'pEntry', 'number', 1);
    end
  end
end

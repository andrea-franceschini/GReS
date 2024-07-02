classdef MeshGlue < handle
  % General mesh gluing class
  % This class reads input file to compute intergrid matrices
  % Mortar method is used to compute the intergrid operator
  % Interfaces structure stores all informations related to each
  % interface

  properties
    interfaces
    model
  end

  methods
    function obj = MeshGlue(modelStruct,fName)
      obj.model = modelStruct;
      readInterfaceFile(obj,fName); % set up all interfaces
    end

    function readInterfaceFile(obj,fName)
      % Set up intergrid operators based on input file
      fID = openReadOnlyFile(fName);
      intStruct = [];
      l = 0;
      while ~feof(fID)
        line = fgetl(fID);
        % read interface datas
        if ~isempty(line) && ~strcmp(line(1),'%')
          str = split(line);
          [d,s] = obj.dealInterfaceLine(str(1:2));
          % create instance of mortar class to compute intergrid operator
          mortar = Mortar3D(1,obj.model(d(1)).Grid.topology,s(1),...
            obj.model(d(2)).Grid.topology,s(2));
          type = str{3};
          switch type
            case 'EB'
              nGP = obj.dealIntegrationParams(fName,d(2),str{4});
              [D,M] = mortar.computeMortarElementBased(nGP);
            case 'RBF'
              [nGP,nInt] = obj.dealIntegrationParams(fName,d(2),str{4:5});
              [D,M] = mortar.computeMortarRBF(nGP,nInt,'gauss');
            otherwise
              error('Invalid tag for integration scheme in %s',fName);
          end
          [~,n_a] = mortar.computeNodalNormal(getElem(mortar,Gauss(mortar.slaveCellType,3,2),'slave'));
          % set interface structure
          intStruct = [intStruct; struct('Master',d(1),'Slave',d(2),...
            'masterSet',mortar.nodesMaster,'slaveSet',mortar.nodesSlave,...
            'nodeNormal',n_a,'InterpOperator',D\M,'mortar',mortar)];
        end
      end
      obj.interfaces = intStruct;
    end

    function [nGP,varargout] = dealIntegrationParams(obj,fName,slaveID,varargin)
      if strcmpi(varargin{1},'default')
        nGP = obj.model(slaveID).Gauss.nNode1D;
      else
        nGP = str2double(varargin{1});
        assert(~isnan(nGP),'Invalid entry for gaussPoints in %s',fName);
      end
      %
      if nargin > 4
        varargout{1} = str2double(varargin{2});
        assert(~isnan(varargout{1}),'Invalid entry for interpolation points in %s',fName);
      end
    end


  end

  methods (Static)
    function [domID,surfID] = dealInterfaceLine(s)
      n1 = sscanf(s{1},'(%i,%i)');
      n2 = sscanf(s{2},'(%i,%i)');
      domID = [n1(1);n2(1)];
      surfID = [n1(2);n2(2)];
    end


  end

end


classdef MeshGlue < handle
   % General mesh gluing class
   % This class implements algorithms to handle non conforming domains
   % Each domain is a separate simulation framework in the model struct
   % This class uses interface information to merge each model into a
   % unique solution system

   properties
      interfaces           % stores interface data
      model                % model setup
      MD_struct            % multidomain info to retrieve linsyst
      domainConn           % domain connectivity
      countDoF
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
                  'nodeNormal',n_a,'InterpOperator',full(D\M),'mortar',mortar,'slaveMat',D,'crossMat',M)];
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


      function obj = setupMDstruct(obj)
         % structure to set up the mortar linear system
         nDom = numel(obj.model); % number of domains
         obj.MD_struct = struct('type',{},'set',{},'tag',{},'dom',{},'physic',{},'entities',{});
         obj.domainConn = cell(nDom,1); % domain connectivity
         nS = 0; %slave interfaces counter
         k = 0;
         for i = 1:nDom
            phList = obj.model(i).DoFManager.subPhysics;
            nPh = numel(phList);
            % get domains in contact
            map = [extractfield(obj.interfaces(:),'Master')',...
               extractfield(obj.interfaces(:),'Slave')'];
            [r,c] = find(map == i);
            id = c == 1; % mark master surfaces
            c_adj = c;
            c_adj(id) = 2; c_adj(~id) = 1;
            obj.domainConn{i} = unique(map(r,c_adj),'stable');
            %
            for ph = 1:nPh
               k = k+1;
               kk = 0;
               inner_nodes = true(obj.model(i).Grid.topology.nNodes,1);
               for j = 1:length(r)
                  if isMortarField(obj,i,r(j),phList(ph))
                     k = k + 1;
                     kk = kk + 1;
                     % check if a physic require mortar integration
                     % store entities list belonging to each set
                     if id(j) % master surface
                        list = obj.interfaces(r(j)).masterSet;
                        tag = 'master';
                     else
                        nS = nS + 1;
                        list = obj.interfaces(r(j)).slaveSet;
                        tag = 'slave';
                     end
                     obj.MD_struct(k).set = obj.model(i).DoFManager.getDoF(phList(ph),list);
                     obj.MD_struct(k).entities = list;
                     obj.MD_struct(k).type = tag;
                     obj.MD_struct(k).tag = r(j);
                     obj.MD_struct(k).dom = i;
                     obj.MD_struct(k).physic = phList(ph);
                     inner_nodes(list) = 0;
                  end
               end
               % store inner nodes
               obj.MD_struct(k-kk).set = obj.model(i).DoFManager.getDoF(phList(ph),find(inner_nodes));
               obj.MD_struct(k-kk).entities = find(inner_nodes);
               obj.MD_struct(k-kk).type = 'inner';
               obj.MD_struct(k-kk).tag = 0;
               obj.MD_struct(k-kk).dom = i;
               obj.MD_struct(k-kk).physic = phList(ph);
            end
         end

         fields = ~strcmp(extractfield(obj.MD_struct,'type'),'slave');
         nFlds = numel(obj.MD_struct);
         cntDoF = zeros(nFlds,1);
         for ii = 1:nFlds
            cntDoF(ii) = numel(obj.MD_struct(ii).set);
         end
         cntDoF(~fields) = 0;
         cntDoF = cumsum(cntDoF);
         obj.countDoF = [0;cntDoF];
      end

      %       function [mat,rhs] = buildMDsyst(obj)
      %           % loop trough the blocks of the final mortar matrix and build blocks
      %           % 1 by 1
      %           flds = find(~strcmp(extractfield(obj.MD_struct,'type'),'slave'));
      %           N = max(obj.countDoF);
      %           rhs = zeros(N,1);
      %           iVec = [];
      %           jVec = [];
      %           mVec = [];
      %           % loop on the matrix
      %           for i = flds
      %               fprintf('processing field %i \n',i)
      %               domID = obj.MD_struct(i).dom;
      %               switch obj.MD_struct(i).type
      %                   % INNER BLOCK ROW
      %                   case 'inner'
      %                       rhs(obj.getDofs_MD(i)) = obj.getRhsBlock(i);
      %                       for j = flds
      %                           if domID == obj.MD_struct(j).dom % same domain
      %                               blk = obj.getBlock(i,j);
      %                               %                         [jvec,ivec] = meshgrid(obj.getDofs_MD(i),obj.getDofs_MD(j));
      %                           else
      %                               if ismember(obj.MD_struct(j).dom,obj.domainConn{domID}) % other connected domain
      %                                   switch obj.MD_struct(j).type
      %                                       case 'inner' % empty block
      %                                           blk = obj.getNullBlock(i,j);
      %                                       case 'master'
      %                                           % find connected slave surface
      %                                           idS = getSlave(obj,j);
      %                                           if domID ~= obj.MD_struct(idS).dom
      %                                               blk = obj.getNullBlock(i,j);
      %                                           else
      %                                               % get mortar operator
      %                                               E = expandMortarOperator(obj,j);
      %                                               blk = obj.getBlock(i,idS)*E;
      %                                           end
      %                                   end
      %                               else % other non-connected domain
      %                                   blk = obj.getNullBlock(i,j);
      %                               end
      %                           end
      %                           [jvec,ivec] = meshgrid(obj.getDofs_MD(i),obj.getDofs_MD(j));
      %                           iVec = [iVec; ivec(:)];
      %                           jVec = [jVec; jvec(:)];
      %                           mVec = [mVec; blk(:)];
      %                       end
      %                       % INTERFACE MASTER BLOCK ROW
      %                   case 'master'
      %                       idS = getSlave(obj,i);
      %                       E = obj.expandMortarOperator(i);
      %                       rhs(obj.getDofs_MD(i)) = obj.getRhsBlock(i) + E'*obj.getRhsBlock(idS);
      %                       for j = flds
      %                           switch obj.MD_struct(j).type
      %                               case 'inner'
      %                                   if domID == obj.MD_struct(j).dom
      %                                       blk = obj.getBlock(i,j);
      %                                   else
      %                                       if ismember(obj.MD_struct(j).dom,obj.domainConn{domID}) % other connected domain
      %                                           % find connected slave surface
      %                                           idS = getSlave(obj,i);
      %                                           if obj.MD_struct(j).dom ~= obj.MD_struct(idS).dom
      %                                               blk = obj.getNullBlock(i,j);
      %                                           else
      %                                               % get mortar operator
      %                                               E = expandMortarOperator(obj,i);
      %                                               blk = E'*obj.getBlock(idS,j);
      %                                           end
      %                                       else
      %                                           blk = obj.getNullBlock(i,j);
      %                                       end
      %                                   end
      %                               case 'master'
      %                                   if domID == obj.MD_struct(j).dom && obj.MD_struct(i).tag == obj.MD_struct(j).tag
      %                                       idS = getSlave(obj,i);
      %                                       E = expandMortarOperator(obj,i);
      %                                       blk = getBlock(obj,i,i)+...
      %                                           E'*getBlock(obj,idS,idS)*E;
      %
      %                                   else % each master interface is disjoint from each other
      %                                       blk = obj.getNullBlock(i,j);
      %                                   end
      %                           end
      %                           [jvec,ivec] = meshgrid(obj.getDofs_MD(i),obj.getDofs_MD(j));
      %                           iVec = [iVec; ivec(:)];
      %                           jVec = [jVec; jvec(:)];
      %                           mVec = [mVec; blk(:)];
      %                       end
      %               end
      %           end
      %           mat = sparse(iVec,jVec,mVec,N,N);
      %           % mat = sparse
      %       end

      function [mat,rhs] = buildMDsyst(obj)
         % loop trough the blocks of the final mortar matrix and build blocks
         % 1 by 1
         flds = find(~strcmp(extractfield(obj.MD_struct,'type'),'slave'));
         nflds = length(flds);
         mat = cell(nflds);
         rhs = cell(nflds,1);
         ii = 1;
         % loop on the matrix
         for i = flds
            %fprintf('Processing field %i \n',ii)
            domID = obj.MD_struct(i).dom;
            jj = 1;
            switch obj.MD_struct(i).type
               % INNER BLOCK ROW
               case 'inner'
                  rhs{ii} = obj.getRhsBlock(i);
                  for j = flds
                     if domID == obj.MD_struct(j).dom % same domain
                        mat{ii,jj} = obj.getBlock(i,j);
                        %                         [jvec,ivec] = meshgrid(obj.getDofs_MD(i),obj.getDofs_MD(j));
                     else
                        if ismember(obj.MD_struct(j).dom,obj.domainConn{domID}) % other connected domain
                           switch obj.MD_struct(j).type
                              case 'inner' % empty block
                                 mat{ii,jj} = obj.getNullBlock(i,j);
                              case 'master'
                                 % find connected slave surface
                                 idS = getSlave(obj,j);
                                 if domID ~= obj.MD_struct(idS).dom
                                    mat{ii,jj} = obj.getNullBlock(i,j);
                                 else
                                    % get mortar operator
                                    E = expandMortarOperator(obj,j);
                                    mat{ii,jj} = obj.getBlock(i,idS)*E;
                                 end
                           end
                        else % other non-connected domain
                           mat{ii,jj} = obj.getNullBlock(i,j);
                        end
                     end
                     jj = jj+1;
                  end
                  % INTERFACE MASTER BLOCK ROW
               case 'master'
                  idS = getSlave(obj,i);
                  E = obj.expandMortarOperator(i);
                  rhs{ii} = obj.getRhsBlock(i) + E'*obj.getRhsBlock(idS);
                  for j = flds
                     switch obj.MD_struct(j).type
                        case 'inner'
                           if domID == obj.MD_struct(j).dom
                              mat{ii,jj} = obj.getBlock(i,j);
                           else
                              if ismember(obj.MD_struct(j).dom,obj.domainConn{domID}) % other connected domain
                                 % find connected slave surface
                                 idS = getSlave(obj,i);
                                 if obj.MD_struct(j).dom ~= obj.MD_struct(idS).dom
                                    mat{ii,jj} = obj.getNullBlock(i,j);
                                 else
                                    % get mortar operator
                                    E = expandMortarOperator(obj,i);
                                    mat{ii,jj} = E'*obj.getBlock(idS,j);
                                 end
                              else
                                 mat{ii,jj} = obj.getNullBlock(i,j);
                              end
                           end
                        case 'master'
                           if domID == obj.MD_struct(j).dom && obj.MD_struct(i).tag == obj.MD_struct(j).tag
                              idS = getSlave(obj,i);
                              E = expandMortarOperator(obj,i);
                              mat{ii,jj} = getBlock(obj,i,i)+...
                                 E'*getBlock(obj,idS,idS)*E;
                           elseif domID == obj.MD_struct(j).dom
                              %mat{ii,jj} = obj.getNullBlock(i,j); % postulate master interface is disjoint from each other
                              mat{ii,jj} = obj.getBlock(i,j);
                           else
                              idS1= getSlave(obj,i);
                              idS2 = getSlave(obj,j);
                              if obj.MD_struct(idS1).dom == obj.MD_struct(idS2).dom % possibly connected slave surface
                                 E1 = expandMortarOperator(obj,i);
                                 E2 = expandMortarOperator(obj,j);
                                 mat{ii,jj} = E1'*obj.getBlock(idS1,idS2)*E2;
                              else
                                 mat{ii,jj} = obj.getNullBlock(i,j);
                              end
                           end
                     end
                     jj = jj+1;
                  end
            end
            ii = ii+1;
         end
      end

      function [J_MD,rhs_MD] = getMDlinSyst(obj)
         setupMDstruct(obj);
         [cellJ,cellRhs] = buildMDsyst(obj);
         [J_MD,rhs_MD] = obj.cell2mat_MG(cellJ,cellRhs);
      end


      function E_exp = expandMortarOperator(obj,id)
         % return mortar operator expanded according to the number of DoFs
         % per node for a certain physic
         nc = obj.model(obj.MD_struct(id).dom).DoFManager.getNumbComp(obj.MD_struct(id).physic);
         E = obj.interfaces(obj.MD_struct(id).tag).InterpOperator;
         idRow = 1:size(E,1);
         idCol = 1:size(E,2);
         E_exp = zeros(nc*size(E,1),nc*size(E,2));
         for k = nc-1:-1:0
            r = nc*idRow - k;
            c = nc*idCol - k;
            E_exp(r,c) = E;
         end
         E_exp = sparse(E_exp);
      end

      function dofs = getDofs_MD(obj,i)
         dofs = 1+obj.countDoF(i):obj.countDoF(i+1);
         dofs = dofs';
      end



      function out = isMortarField(obj,domID,intID,phID)
         % check if a certain physic requires mortar interpolation for a
         % certain domain
         % get domain in contacts
         interfMap = [extractfield(obj.interfaces(:),'Master')',...
            extractfield(obj.interfaces(:),'Slave')'];
         domConn = interfMap(intID,interfMap(intID,:) ~= domID);
         % domConn = tmp(intID);
         %[r,~] = find(interfMap == domID);
         % get connected domain
         %domConn = interfMap(interfMap(r(intID),:)~=domID);
         % find all connected domains
         out = ismember(phID,obj.model(domConn).DoFManager.subPhysics);
      end

      function blk = getNullBlock(obj,i,j)
         nr = length(obj.MD_struct(i).set);
         nc = length(obj.MD_struct(j).set);
         blk = sparse(nr,nc);
      end

      function blk = getBlock(obj,i,j)
         r = obj.MD_struct(i).set;
         c = obj.MD_struct(j).set;
         domID = obj.MD_struct(i).dom;
         J = cell2mat(obj.model(domID).Discretizer.J);
         blk = J(r,c);
      end

      function blk = getRhsBlock(obj,i)
         r = obj.MD_struct(i).set;
         domID = obj.MD_struct(i).dom;
         rhs = cell2mat(obj.model(domID).Discretizer.rhs);
         blk = rhs(r);
      end

      function idSlave = getSlave(obj,id)
         % get slave set connected to a master interface
         assert(strcmp(obj.MD_struct(id).type,'master'),['Input ' ...
            'interface is not a master set']);
         % get interface tag
         tag = obj.MD_struct(id).tag;
         idSlave = find(all([extractfield(obj.MD_struct,'tag')'==tag,...
            strcmp(extractfield(obj.MD_struct,'type')','slave')],2));
      end

      function idMaster = getMaster(obj,id)
         % get slave set connected to a master interface
         assert(strcmp(obj.MD_struct(id).type,'slave'),['Input ' ...
            'interface is not a slave set']);
         % get interface tag
         tag = obj.MD_struct(id).tag;
         idMaster = find(all([extractfield(obj.MD_struct,'tag')'==tag,...
            strcmp(extractfield(obj.MD_struct,'type')','master')],2));
      end


   end

   methods (Static)
      function [domID,surfID] = dealInterfaceLine(s)
         n1 = sscanf(s{1},'(%i,%i)');
         n2 = sscanf(s{2},'(%i,%i)');
         domID = [n1(1);n2(1)];
         surfID = [n1(2);n2(2)];
      end


      function [J_mat,rhs_mat] = cell2mat_MG(cellJ,cellRhs)
         J_mat = [];
         N = size(cellJ,1);
         rhs_mat = cell2mat(cellRhs);
         clear cellRhs
         for i = 1:N
            J_mat = [J_mat; cell2mat(cellJ(i,:))];
            cellJ(i,:) = cell(1,N);
         end
      end


   end

end


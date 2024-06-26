function modStruct = buildModelStruct(fName)
% build a structure array containing istances of all classes for different
% domains
fID = openReadOnlyFile(fName);
dofFile = []; bcList = []; bc = [];
l = 0; c = 0;
modStruct = [];
while ~feof(fID)
   line = readToken(fID,fName);
   l = l+1;
   if ~isempty(line) && ~strcmp(line(1),'%')
      c = c+1;
      nextline = readToken(fID,fName);
      l = l+1;
      while ~strcmp(nextline,'<EndDomain>')
         switch nextline
            case '<Name>'
               name = readToken(fID,fName);
               l = l+1;
            case '<ModelType>'
               mods = convertCharsToStrings(readToken(fID,fName));
               l = l+1;
               mods = (split(mods))';
            case '<SimulationParameters>'
               sim = readToken(fID,fName);
               l = l+1;
            case '<Geometry>'
               topology = Mesh();
               mshline = readToken(fID,fName);
               l = l+1;
               if contains(mshline,'.msh')
                  topology.importGMSHmesh(mshline);
               elseif contains(mshline,'.vtk')
                  topology.importVTKmesh(Mesh(),mshline);
               end
            case '<Gauss>'
               p = fscanf(fID,'%f');
               gauss = Gauss(p(1),p(2),p(3));
            case '<Material>'
               matFile = readToken(fID,fName);
               l = l+1;
            case '<DoFManager>'
               dofFile = readToken(fID,fName);
               l = l+1;
            case '<OutState>'
               outFile = readToken(fID,fName);
               l = l+1;
            case '<BoundaryConditions>'
               nBC = fscanf(fID,'%i');
               assert(~isempty(nBC),['Expected number of BCs as',...
                  'first row of <BoundaryConditions> block in %s file'],fName);
               l = l+1;
               bcList = strings(1,nBC);
               for i = 1:nBC
                  bcName = readToken(fID,fName);
                  assert(~strcmp(bcName([1,end]),'<>'),['Incorrect number of BCs' ...
                     ' in <BoundaryConditions block of %s file'],fName);
                  bcList(i) = bcName;
                  l = l+1;
               end
            otherwise
               if ~isempty(nextline)
                  error('Invalid block name in %s domain input file',fName);
               else
                  error(['Error in %s input file: Invalid blank line, ...' ...
                     'in <Domain> block'],fName);
               end
         end
         nextline = readToken(fID,fName);
         l = l+1;
         assert(strcmp(nextline([1,end]),'<>'),'Syntax error or unexpected text in line %i of file %s',l,fName);
      end
      model = ModelType(mods);
      simParam = SimulationParameters(model,sim);
      mat = Materials(model,matFile);
      if ~isempty(dofFile)
         dof = DoFManager(topology,model,dofFile);
      else
         dof = DoFManager(topology,model);
      end
      elems = Elements(topology,gauss);
      faces = Faces(model,topology);
      grid = struct('topology',topology,'cells',elems,'faces',faces);
      if isempty(bcList)
         warning('Undefined Boundary conditions for Domain %i in file %s',c,fName);
      else
         bc = Boundaries(bcList,model,grid,dof);
      end
      state = State(model,grid,mat,gauss);
      printUtils = OutState(model,mat,grid,outFile,name,gauss);
      linSyst = Discretizer(model,simParam,dof,grid,mat,gauss);
      modStruct = [modStruct;struct('id',c,'DomainName',name,'ModelType',model,...
         'SimParams',simParam,'Grid',grid,'Material',mat,'Gauss',gauss,...
         'DoFManager',dof,'BoundaryConditions',bc,'State',state,'OutState',printUtils,...
         'Discretizer',linSyst)];
   end
end
end


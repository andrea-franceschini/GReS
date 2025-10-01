function writeBCfiles(fName,item,type,physic,bcName,time,vals,varargin)
% utility to write BC input files (only constant BC for each time step are
% allowed) Entities are identified by surface ID or by direct assignment
% Example
% writeBCfiles('BCs/dirFlow','SurfBC','Dir','Flow','flowBCname',0,0,mesh,1)
% writeBCfiles('BCs/dirPoroXY','SurfBC','Dir',{Poro,x,y},'poroBCname',0,0,mesh,1)

assert(length(time)==length(vals),'BC times and values sets must have equal size');

% Writing general BCs file
fID = fopen(strcat(fName,'.dat'),'w');

% item ('NodeBC','SurfBC','VolumeForce')
fprintf(fID,'%s            %% BC item \n',item);

if ~strcmp(item,'VolumeForce')
  fprintf(fID,'%s            %% BC type \n',type);
end

% Physic
physic = string(physic);
if numel(physic)>1
  assert(strcmp(physic(1),'Poromechanics'),['Direction specification is allowed' ...
    'only for Poromechanics']);
  ph = physic(1);
else
  ph = physic;
end

fprintf(fID,'%s            %% Physics \n',physic(1));

% Direction
if strcmp(physic(1),'Poromechanics')
  dir = physic(2:end);
  %    if strcmp(type,'Neu') && strcmp(item,'SurfBC')
  %       assert(numel(physic)==2,['Only one direction at time is allowed for' ...
  %          ' Poromechanics Neumann BCs ']);
  %       fprintf(fID,'%s \n',dir);
  %    end
end

% BC Name
fprintf(fID,'%s            %% BC name \n',bcName);

% BC list file name
listName = strcat(fName,'/list');
fprintf(fID,'%s \n',listName);

% BC time file
for i = 0:length(time)-1
  fprintf(fID,'%2.6f %s/time%i.dat \n',time(i+1),fName,i);
end

% End file
fprintf(fID,'End');

if ~isfolder(fName)
  mkdir(fName);
end

% Writing BC list of constrained entities
assert(numel(varargin)>0,'Not enough input arguments')
fList = fopen(listName,'w');
if ~isa(varargin{1},"Mesh") % direct assignment
  list = varargin{1};
  if numel(varargin)>1
    vals = varargin{2}; % overwrite standard vals
    assert(length(list)==length(vals))
  end
else
  surf = find(ismember(varargin{1}.surfaceTag,varargin{2}));
  switch item
    case 'NodeBC'
      nodes = unique(varargin{1}.surfaces(surf,:));
      list = nodes;
    case 'SurfBC'
      list = surf;
  end
end

if ismember(ph,["SinglePhaseFlow","Poisson","VariablySaturatedFlow"]) || strcmp(item,'VolumeForce')
  fprintf(fList,'%i         %% Number of fixed entities \n',length(list));
  fprintf(fList,'%i \n',list);
else
  tmp = ismember(["x","y","z"],dir);
  fprintf(fList,'%i ',tmp*length(list));
  fprintf(fList,'   %% Number of fixed entities \n');
  list = repmat(list,sum(tmp),1);
  fprintf(fList,'%i \n',list);
end

% Writing BC vals for each time step
for i = 1:length(time)
  t_name = strcat(fName,'/time',num2str(i-1),'.dat');
  ft = fopen(t_name,'w');
  fprintf(ft,'%%Time %2.4f \n',time(i));
  if length(vals) == length(time)
    if (isscalar(vals(i)))
    fprintf(ft,'%1.6e \n',repelem(vals(i),length(list)));
    end
  else
    fprintf(ft,'%1.6e \n',vals);
  end
end
end
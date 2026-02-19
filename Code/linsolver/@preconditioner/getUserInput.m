function [params] = getUserInput(obj,params,usrInput)


if isempty(usrInput)
  gresLog().log(3,'Using default values for preconditioner\n');
  return
end

usrInput = readstruct(usrInput,AttributeSuffix="");
input = [];

% Unpack the fields and select the correct one
if obj.phys == 0
  if isfield(usrInput,'Flow')
    input = usrInput.Flow;
  end
else
  if isfield(usrInput,'Mechanics')
    input = usrInput.Mechanics;
  end
end

% If the file does not exist use the defaults
if ~isempty(input)
  if obj.DEBUGflag
    fprintf('Using user defined values for preconditioner\n');
  end

  % Get amg params
  if isfield(input,'amg')
    amg = input.('amg');
    params.amg.nLevMax     = getXMLData(amg,params.amg.nLevMax,'nLevMax');
    params.amg.maxCoarseSZ = getXMLData(amg,params.amg.maxCoarseSZ,'maxCoarseSZ');
  end

  % Get smoother params
  if isfield(input,'smoother')
    smoother = input.('smoother');
    params.smoother.nstep     = getXMLData(smoother,params.smoother.nstep,'nstep');
    params.smoother.epsilon   = getXMLData(smoother,params.smoother.epsilon,'epsilon');
    params.smoother.step_size = getXMLData(smoother,params.smoother.step_size,'step_size');
  end

  % Get prolong params
  if isfield(input,'prolong')
    prolong = input.('prolong');
    params.prolong.prol_emin = getXMLData(prolong,params.prolong.prol_emin,'prol_emin');
  end

  % Get coarsen params
  if isfield(input,'coarsen')
    coarsen = input.('coarsen');
    params.coarsen.tau      = getXMLData(coarsen,params.coarsen.tau,'tau');
    params.coarsen.nl_agg   = getXMLData(coarsen,params.coarsen.nl_agg,'nl_agg');
    params.coarsen.SoC_type = getXMLData(coarsen,params.coarsen.SoC_type,'SoC_type');
  end

  % Get test space params
  if isfield(input,'tspace')
    tspace = input.('tspace');
    params.tspace.ntv    = getXMLData(tspace,params.tspace.ntv,'ntv');
    params.tspace.itmax  = getXMLData(tspace,params.tspace.itmax,'itmax');
    params.tspace.method = getXMLData(tspace,params.tspace.method,'method');
  end

else
  if obj.DEBUGflag
    fprintf('Using default values for preconditioner\n');
  end
end
end


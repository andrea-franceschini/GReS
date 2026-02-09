function [params] = getUserInput(obj,params,input)


% If the file does not exist use the defaults
if ~isempty(input)
  if gresLog().getVerbosity > 2
    fprintf('Using user defined values for preconditioner\n');
  end


  if isfield(input,'amg')
    amg = input.('amg');
    params.amg.nLevMax     = getinputData(amg,params.amg.nLevMax,'nLevMax');
    params.amg.maxCoarseSZ = getinputData(amg,params.amg.maxCoarseSZ,'maxCoarseSZ');
  end

  if isfield(input,'smoother')
    smoother = input.('smoother');
    params.smoother.nstep     = getinputData(smoother,params.smoother.nstep,'nstep');
    params.smoother.epsilon   = getinputData(smoother,params.smoother.epsilon,'epsilon');
    params.smoother.step_size = getinputData(smoother,params.smoother.step_size,'step_size');
  end

  if isfield(input,'prolong')
    prolong = input.('prolong');
    params.prolong.prol_emin = getinputData(prolong,params.prolong.prol_emin,'prol_emin');
  end

  if isfield(input,'coarsen')
    coarsen = input.('coarsen');
    params.coarsen.tau      = getinputData(coarsen,params.coarsen.tau,'tau');
    params.coarsen.nl_agg   = getinputData(coarsen,params.coarsen.nl_agg,'nl_agg');
    params.coarsen.SoC_type = getinputData(coarsen,params.coarsen.SoC_type,'SoC_type');
  end

  if isfield(input,'tspace')
    tspace = input.('tspace');
    params.tspace.ntv    = getinputData(tspace,params.tspace.ntv,'ntv');
    params.tspace.itmax  = getinputData(tspace,params.tspace.itmax,'itmax');
    params.tspace.method = getinputData(tspace,params.tspace.method,'method');
  end

else
  if gresLog().getVerbosity > 2
    fprintf('Using default values for preconditioner\n');
  end
end
end


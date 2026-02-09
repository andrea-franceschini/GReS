function [params] = getUserInput(obj,params,xml)
   
   % If the file does not exist use the defaults
   if isfile(xml)
      if gresLog().getVerbosity > 2
         fprintf('Using user defined values for preconditioner\n');
      end

      input = readstruct(xml,AttributeSuffix="");

      if isfield(input,'amg')
         amg = input.('amg');
         params.amg.nLevMax     = getXMLData(amg,params.amg.nLevMax,'nLevMax');
         params.amg.maxCoarseSZ = getXMLData(amg,params.amg.maxCoarseSZ,'maxCoarseSZ');
      end

      if isfield(input,'smoother')
         smoother = input.('smoother');
         params.smoother.nstep     = getXMLData(smoother,params.smoother.nstep,'nstep');
         params.smoother.epsilon   = getXMLData(smoother,params.smoother.epsilon,'epsilon');
         params.smoother.step_size = getXMLData(smoother,params.smoother.step_size,'step_size');
      end
      
      if isfield(input,'prolong')
         prolong = input.('prolong');
         params.prolong.prol_emin = getXMLData(prolong,params.prolong.prol_emin,'prol_emin');
      end

      if isfield(input,'coarsen')
         coarsen = input.('coarsen');
         params.coarsen.tau      = getXMLData(coarsen,params.coarsen.tau,'tau');
         params.coarsen.nl_agg   = getXMLData(coarsen,params.coarsen.nl_agg,'nl_agg');
         params.coarsen.SoC_type = getXMLData(coarsen,params.coarsen.SoC_type,'SoC_type');
      end

      if isfield(input,'tspace')
         tspace = input.('tspace');
         params.tspace.ntv    = getXMLData(tspace,params.tspace.ntv,'ntv');
         params.tspace.itmax  = getXMLData(tspace,params.tspace.itmax,'itmax');
         params.tspace.method = getXMLData(tspace,params.tspace.method,'method');
      end

   else
      if gresLog().getVerbosity > 2
         fprintf('Using default values for preconditioner\n');
      end
   end
end


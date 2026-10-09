function [params] = getUserInput(obj,params,usrInput)

   input = [];

   % Unpack the fields and select the correct one
   if ~isempty(obj.phys) && obj.phys == 0
      if isfield(usrInput,'Flow')
         input = usrInput.Flow;
      end
   elseif ~isempty(obj.phys) && obj.phys >= 1
      if isfield(usrInput,'Mechanics')
         input = usrInput.Mechanics;
      end
   else
      if isfield(usrInput,'Mechanics')
         input = usrInput.Mechanics;
      elseif isfield(usrInput,'Flow')
         input = usrInput.Flow;
      end
   end
   
   % If user parameters are provided, merge allowed parameters
   if ~isempty(input)
      gresLog().log(3,'Using user defined values for preconditioner\n');

      % Lists of parameters blocked from XML configuration
      blockedSmoother = {'nupre', 'nupost', 'nthread', 'method'};
      blockedProlong = {'np', 'itmax_Vol', 'tol_Vol', 'dist_min', 'dist_max', ...
                        'maxcond', 'maxrownrm', 'eps_prol', 'updateCF', ...
                        'patt_pow', 'patt_tau', 'nnzr_max', 'itmax_emin', ...
                        'condmax_emin', 'prec_emin', 'solv_emin', ...
                        'min_lfil', 'max_lfil', 'D_lfil'};
      blockedTspace = {'init_approx'};
   
      % Get amg params
      if isfield(input,'amg')
         params.amg = readInput(params.amg,input.amg);
      end
   
      % Get smoother params (with blocked fields stripped)
      if isfield(input,'smoother')
         smootherInput = removeBlocked(input.smoother, blockedSmoother);
         params.smoother = readInput(params.smoother,smootherInput);
      end
   
      % Get prolong params (with blocked fields stripped)
      if isfield(input,'prolong')
         prolongInput = removeBlocked(input.prolong, blockedProlong);
         params.prolong = readInput(params.prolong,prolongInput);
      end
   
      % Get coarsen params
      if isfield(input,'coarsen')
         params.coarsen = readInput(params.coarsen,input.coarsen);
      end
   
      % Get test space params (with blocked fields stripped)
      if isfield(input,'tspace')
         tspaceInput = removeBlocked(input.tspace, blockedTspace);
         params.tspace = readInput(params.tspace,tspaceInput);
      end
   
   else
      gresLog().log(3,'Using default values for preconditioner\n');
   end
end

function s = removeBlocked(s, blockedList)
   if ~isstruct(s)
      return;
   end
   fnames = fieldnames(s);
   for i = 1:numel(fnames)
      if any(strcmpi(fnames{i}, blockedList))
         s = rmfield(s, fnames{i});
      end
   end
end


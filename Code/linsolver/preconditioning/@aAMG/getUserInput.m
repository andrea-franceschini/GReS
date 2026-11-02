function [params] = getUserInput(obj,params,usrInput)

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
      gresLog().log(3,'Using user defined values for preconditioner\n');
   
      % Get amg params
      if isfield(input,'amg')
         params.amg = readInput(params.amg,input.amg);
      end
   
      % Get smoother params
      if isfield(input,'smoother')
         params.smoother = readInput(params.smoother,input.smoother);
      end
   
      % Get prolong params
      if isfield(input,'prolong')
         params.prolong = readInput(params.prolong,input.prolong);
      end
   
      % Get coarsen params
      if isfield(input,'coarsen')
         params.coarsen = readInput(params.coarsen,input.coarsen);
      end
   
      % Get test space params
      if isfield(input,'tspace')
         params.tspace = readInput(params.tspace,input.tspace);
      end
   
   else
      gresLog().log(3,'Using default values for preconditioner\n');
   end
end


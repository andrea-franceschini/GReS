function [param] = getUserInput(obj,param,usrInput)

   input = [];

   % Unpack the fields and select the correct one based on physics
   if isfield(usrInput, 'fsai')
      input = usrInput.fsai;
   elseif ~isempty(obj.phys) && obj.phys == 0
      if isfield(usrInput, 'Flow') && isfield(usrInput.Flow, 'fsai')
         input = usrInput.Flow.fsai;
      end
   elseif ~isempty(obj.phys) && obj.phys >= 1
      if isfield(usrInput, 'Mechanics') && isfield(usrInput.Mechanics, 'fsai')
         input = usrInput.Mechanics.fsai;
      end
   end

   % Fallback if input not found under Flow/Mechanics
   if isempty(input)
      if isfield(usrInput, 'fsai')
         input = usrInput.fsai;
      end
   end

   if ~isempty(input)
      gresLog().log(3, 'Using user defined values for FSAI preconditioner\n');
      if isfield(input, 'nstep')
         param.nstep = str2double(string(input.nstep));
      end
      if isfield(input, 'step_size')
         param.step_size = str2double(string(input.step_size));
      end
      if isfield(input, 'epsilon')
         param.epsilon = str2double(string(input.epsilon));
      end
      % Note: 'method' is chosen automatically based on symmetry and cannot be set from XML
   else
      gresLog().log(3, 'Using default values for FSAI preconditioner\n');
   end
end

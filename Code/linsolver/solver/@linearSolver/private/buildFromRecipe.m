function prec = buildFromRecipe(obj, recipe, generalsolver)
% BUILDFROMRECIPE Recursively instantiates preconditioners from declarative recipe struct
%
% Usage:
%   prec = obj.buildFromRecipe(recipe, generalsolver)

   if ischar(recipe) || isstring(recipe)
      recipe = struct('type', string(recipe));
   end

   debugflag = obj.DEBUGflag;

   switch string(recipe.type)
      case "aAMG"
         phys = getFieldOr(recipe, 'physics', "displacements");
         prec = aAMG(debugflag, generalsolver, phys);
      case "multiDomAMG"
         phys = getFieldOr(recipe, 'physics', "pressure");
         prec = multiDomAMG(debugflag, generalsolver, phys);
      case "RACP"
         phys = getFieldOr(recipe, 'physics', "displacements");
         prec = RACP(debugflag, generalsolver, phys);
      case "fixedStress"
         prec = fixedStress(debugflag, generalsolver);
      case "CPR"
         prec = CPR(debugflag, generalsolver);
      case "efemPrec"
         prec = efemPrec(debugflag, generalsolver, obj.nsyTol);
      otherwise
         gresLog().warning(3, sprintf('Unknown preconditioner type: %s. Falling back to MATLAB solver.', recipe.type));
         prec = [];
         return;
   end

   % Recursively build and inject child subpreconditioners
   if isfield(recipe, 'subPrecs') && isstruct(recipe.subPrecs)
      fields = fieldnames(recipe.subPrecs);
      args = cell(1, 2 * numel(fields));
      for i = 1:numel(fields)
         fName = fields{i};
         childNode = recipe.subPrecs.(fName);
         childPrec = buildFromRecipe(obj, childNode, generalsolver);
         if isempty(childPrec)
            prec = [];
            return;
         end
         args{2*i - 1} = fName;
         args{2*i}     = childPrec;
      end
      prec.setSubPrecs(args{:});
   end
end


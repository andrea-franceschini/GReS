function prec = resolvePreconditioner(obj, generalsolver, activeBlocks, order)
% RESOLVEPRECONDITIONER Recursively constructs composite preconditioners for active blocks
%
% Usage:
%   prec = obj.resolvePreconditioner(generalsolver, activeBlocks, order)
%
% Supported orders:
%   - "physicsFirst" (default): physics coupling outer, domain/interfaces inner
%   - "racpFirst": saddle-point interface outer (RACP), physics coupling inner

   if nargin < 4 || isempty(order)
      order = "physicsFirst";
   end
   if nargin < 3 || isempty(activeBlocks)
      activeBlocks = 1:obj.systemInfo.nBlocks;
   end

   blocks = obj.systemInfo.blocks(activeBlocks);
   uniqueVars = unique([blocks.varName]);
   debugflag = obj.DEBUGflag;

   % -------------------------------------------------------------
   % BASE CASE: Single physical variable remaining in active scope
   % -------------------------------------------------------------
   if numel(uniqueVars) == 1
      [chronosFlag, prec] = chooseSinglePhys(obj, generalsolver, debugflag, uniqueVars, activeBlocks);
      if ~chronosFlag
         prec = [];
      end
      return;
   end

   % -------------------------------------------------------------
   % TIER: Saddle-Point Lagrange Multipliers
   % -------------------------------------------------------------
   isSaddle = [blocks.kind] == "interface" & ...
              contains([blocks.solverType], ["MeshTying", "Contact"]);
   if strcmpi(order, "racpFirst") && any(isSaddle)
      % Isolate physical variables and non-saddle interfaces into A11
      varBlockIds = activeBlocks(~isSaddle);

      innerPrec = resolvePreconditioner(obj, generalsolver, varBlockIds, order);
      if isempty(innerPrec)
         prec = [];
         return;
      end

      prec = RACP(debugflag, generalsolver, uniqueVars(1));
      prec.setSubPrecs(innerPrec);
      return;
   end

   % -------------------------------------------------------------
   % TIER: Poromechanics Biot Coupling (displacements + flow)
   % -------------------------------------------------------------
   if ismember("displacements", uniqueVars) && any(ismember(["pressure", "saturation"], uniqueVars))
      uBlockIds = activeBlocks(ismember([blocks.varName], "displacements"));
      flowBlockIds = activeBlocks(ismember([blocks.varName], ["pressure", "saturation"]));

      mechPrec = resolvePreconditioner(obj, generalsolver, uBlockIds, order);
      flowPrec = resolvePreconditioner(obj, generalsolver, flowBlockIds, order);

      if isempty(mechPrec) || isempty(flowPrec)
         prec = [];
         return;
      end

      prec = fixedStress(debugflag, generalsolver);
      prec.setSubPrecs("mechanics", mechPrec, "flow", flowPrec);
      return;
   end

   % -------------------------------------------------------------
   % TIER: Two-Phase Flow Coupling (pressure + saturation)
   % -------------------------------------------------------------
   if ismember("pressure", uniqueVars) && ismember("saturation", uniqueVars)
      pBlockIds = activeBlocks(ismember([blocks.varName], "pressure"));
      pPrec = resolvePreconditioner(obj, generalsolver, pBlockIds, order);

      if isempty(pPrec)
         prec = [];
         return;
      end

      prec = CPR(debugflag, generalsolver);
      prec.setSubPrecs("pressure", pPrec);
      return;
   end

   % -------------------------------------------------------------
   % TIER: Embedded Fracture Coupling (displacements + fractureJump)
   % -------------------------------------------------------------
   if ismember("displacements", uniqueVars) && ismember("fractureJump", uniqueVars)
      uBlockIds = activeBlocks(ismember([blocks.varName], "displacements"));
      mechPrec = resolvePreconditioner(obj, generalsolver, uBlockIds, order);

      if isempty(mechPrec)
         prec = [];
         return;
      end

      prec = efemPrec(debugflag, generalsolver, obj.nsyTol);
      prec.setSubPrecs("displacements", mechPrec);
      return;
   end

   % Fallback: could not resolve
   gresLog().warning(3, sprintf('No preconditioner combination available for active fields: %s. Falling back to MATLAB solver.', strjoin(uniqueVars, ', ')));
   prec = [];
end

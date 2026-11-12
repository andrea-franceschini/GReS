function info = describeSystem(generalsolver, verbose)
% DESCRIBESYSTEM Inspect generalsolver and describe the assembled matrix structure
%
% Usage:
%   info = linearSolver.describeSystem(generalsolver, verbose)
%
% Description:
%   Inspects the domains, variables, and interfaces configured in generalsolver
%   and details the global block matrix layout (domains followed by interface
%   multipliers), including intra-domain physics couplings (e.g., Biot coupling)
%   and inter-domain couplings (e.g., mesh tying, contact).

   if nargin < 2
      verbose = true;
   end

   nDom = generalsolver.nDom;
   nInterf = generalsolver.nInterf;
   domains = generalsolver.domains;
   interfaces = generalsolver.interfaces;

   % Count variables and collect metadata per domain variable block
   varBlocks = struct('blockIdx', {}, 'domainId', {}, 'varId', {}, ...
                      'varName', {}, 'nDof', {}, 'solvers', {});
   domVarOffsets = zeros(nDom, 1);
   totalVars = 0;

   for iD = 1:nDom
      dom = domains(iD);
      dofm = dom.dofm;
      nV = dofm.getNumberOfVariables();
      varNames = dofm.getVariableNames();
      domVarOffsets(iD) = totalVars;

      for v = 1:nV
         varName = varNames(v);
         nDof = dofm.getNumbDoF(varName);

         % Detect which physics solvers are associated with this variable
         solvs = [];
         for sIdx = 1:numel(dom.solverNames)
            sName = dom.solverNames(sIdx);
            solvObj = dom.getPhysicsSolver(sName);
            if any(strcmpi(solvObj.getField(), varName))
               solvs = [solvs, sName]; %#ok<AGROW>
            end
         end

         blockIdx = totalVars + v;
         varBlocks(blockIdx).blockIdx = blockIdx;
         varBlocks(blockIdx).domainId = iD;
         varBlocks(blockIdx).varId = v;
         varBlocks(blockIdx).varName = string(varName);
         varBlocks(blockIdx).nDof = nDof;
         varBlocks(blockIdx).solvers = solvs;
      end
      totalVars = totalVars + nV;
   end

   % Collect metadata per interface block
   interfBlocks = struct('blockIdx', {}, 'interfId', {}, 'interfType', {}, ...
                         'multiplierLocation', {}, 'nDof', {}, ...
                         'coupledVar', {}, 'slaveDom', {}, 'masterDom', {}, ...
                         'slaveBlock', {}, 'masterBlock', {});

   for iI = 1:nInterf
      interf = interfaces{iI};
      blockIdx = totalVars + iI;
      sDomId = interf.domainId(1);
      mDomId = interf.domainId(2);
      cVar = interf.getCoupledVariables();

      sDofm = domains(sDomId).dofm;
      mDofm = domains(mDomId).dofm;
      sVarId = sDofm.getVariableId(cVar);
      mVarId = mDofm.getVariableId(cVar);

      sBlock = domVarOffsets(sDomId) + sVarId;
      mBlock = domVarOffsets(mDomId) + mVarId;

      interfBlocks(iI).blockIdx = blockIdx;
      interfBlocks(iI).interfId = iI;
      interfBlocks(iI).interfType = class(interf);
      interfBlocks(iI).multiplierLocation = char(interf.multiplierLocation);
      interfBlocks(iI).nDof = interf.getNumbDoF();
      interfBlocks(iI).coupledVar = string(cVar);
      interfBlocks(iI).slaveDom = sDomId;
      interfBlocks(iI).masterDom = mDomId;
      interfBlocks(iI).slaveBlock = sBlock;
      interfBlocks(iI).masterBlock = mBlock;
   end

   % Unified diagonal descriptor array
   nTotalBlocks = totalVars + nInterf;
   blocks = struct('blockIdx', {}, 'kind', {}, 'domainId', {}, ...
                   'varName', {}, 'nDof', {}, 'solverType', {});
   for b = 1:totalVars
      vb = varBlocks(b);
      blocks(b).blockIdx = b;
      blocks(b).kind = "variable";
      blocks(b).domainId = vb.domainId;
      blocks(b).varName = vb.varName;
      blocks(b).nDof = vb.nDof;
      blocks(b).solverType = vb.solvers;
   end
   for iI = 1:nInterf
      ib = interfBlocks(iI);
      b = totalVars + iI;
      blocks(b).blockIdx = b;
      blocks(b).kind = "interface";
      blocks(b).domainId = [ib.slaveDom, ib.masterDom];
      blocks(b).varName = ib.coupledVar;
      blocks(b).nDof = ib.nDof;
      blocks(b).solverType = ib.interfType;
   end

   % Identify couplings between blocks
   couplings = {};

   % 1. Intra-domain couplings (e.g. Biot, MultiphasePoromechanics, EFEM)
   for iD = 1:nDom
      dom = domains(iD);
      dofm = dom.dofm;
      offset = domVarOffsets(iD);

      for sIdx = 1:numel(dom.solverNames)
         sName = dom.solverNames(sIdx);
         solvObj = dom.getPhysicsSolver(sName);

         if isa(solvObj, 'BiotFullyCoupled')
            uBlock = offset + dofm.getVariableId('displacements');
            pBlock = offset + dofm.getVariableId('pressure');
            couplings{end+1} = sprintf('Blocks (%d,%d) and (%d,%d): Domain %d coupled with Biot (%s)', ...
               uBlock, uBlock, pBlock, pBlock, iD, sName); %#ok<AGROW>

            if isa(solvObj, 'MultiphasePoromechanics') && dofm.isVariable('saturation')
               sBlock = offset + dofm.getVariableId('saturation');
               couplings{end+1} = sprintf('Blocks (%d,%d) and (%d,%d): Domain %d coupled with multiphase flow/saturation', ...
                  pBlock, pBlock, sBlock, sBlock, iD); %#ok<AGROW>
            end
         elseif isa(solvObj, 'EmbeddedFractureMechanics') && dofm.isVariable('fractureJump')
            uBlock = offset + dofm.getVariableId('displacements');
            jBlock = offset + dofm.getVariableId('fractureJump');
            couplings{end+1} = sprintf('Blocks (%d,%d) and (%d,%d): Domain %d embedded fracture mechanics coupling', ...
               uBlock, uBlock, jBlock, jBlock, iD); %#ok<AGROW>
         end
      end
   end

   % 2. Inter-domain couplings via interfaces
   for iI = 1:nInterf
      ib = interfBlocks(iI);
      desc = ib.interfType;
      if contains(desc, 'Contact')
         descType = 'contact interface';
      elseif contains(desc, 'MeshTying')
         descType = 'mesh tying';
      else
         descType = sprintf('interface (%s)', desc);
      end

      couplings{end+1} = sprintf('Block (%d,%d) [%s]: couples block (%d,%d) (Domain %d, %s) with block (%d,%d) (Domain %d, %s) via %s', ...
         ib.blockIdx, ib.blockIdx, desc, ...
         ib.slaveBlock, ib.slaveBlock, ib.slaveDom, ib.coupledVar, ...
         ib.masterBlock, ib.masterBlock, ib.masterDom, ib.coupledVar, ...
         descType); %#ok<AGROW>
   end

   % Build structured output
   info.nBlocks = nTotalBlocks;
   info.nDomainVariables = totalVars;
   info.nInterfaces = nInterf;
   info.blocks = blocks;
   info.varBlocks = varBlocks;
   info.interfBlocks = interfBlocks;
   info.couplings = couplings;

   % Print summary to console when requested
   if verbose
      fprintf('=== Assembled Matrix Structure (%dx%d blocks) ===\n', nTotalBlocks, nTotalBlocks);
      for b = 1:totalVars
         vb = varBlocks(b);
         solvStr = join(vb.solvers, ',');
         fprintf('Block (%d,%d): Domain %d, variable "%s" [DoFs: %d, physics: %s]\n', ...
            b, b, vb.domainId, vb.varName, vb.nDof, solvStr);
      end
      for iI = 1:nInterf
         ib = interfBlocks(iI);
         fprintf('Block (%d,%d): Interface %d (%s) [DoFs: %d, location: %s, variable: "%s"]\n', ...
            ib.blockIdx, ib.blockIdx, ib.interfId, ib.interfType, ib.nDof, ib.multiplierLocation, ib.coupledVar);
      end

      if ~isempty(couplings)
         fprintf('\nCouplings:\n');
         for c = 1:numel(couplings)
            fprintf(' - %s\n', couplings{c});
         end
      end
      fprintf('===================================================\n');
   end
end

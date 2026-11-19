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

   % Loop over the domains
   for iD = 1:nDom
      dom = domains(iD);
      dofm = dom.dofm;
      nV = dofm.getNumberOfVariables();
      varNames = dofm.getVariableNames();
      domVarOffsets(iD) = totalVars;

      % Loop over the variables in the domain iD
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

   % Loop over all the interfaces
   for iI = 1:nInterf
      interf = interfaces{iI};
      blockIdx = totalVars + iI;
      
      % Get master and slave domains
      sDomId = interf.domainId(1);
      mDomId = interf.domainId(2);
      cVar = interf.getCoupledVariables();

      % Get master and slave domains degrees of freedom
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
                   'varName', {}, 'nDof', {}, 'solverType', {}, 'interfId',{});
   
   % Loop over all the variables
   for b = 1:totalVars
      vb = varBlocks(b);
      blocks(b).blockIdx = b;
      blocks(b).kind = "variable";
      blocks(b).domainId = vb.domainId;
      blocks(b).varName = vb.varName;
      blocks(b).nDof = vb.nDof;
      blocks(b).solverType = vb.solvers;
   end
   
   % Loop over all the interface blocks
   for iI = 1:nInterf
      ib = interfBlocks(iI);
      b = totalVars + iI;
      blocks(b).blockIdx = b;
      
      % At the moment MeshTyingTPFA is a particular interface as it couples
      % pressures with pressures through other pressures so treated as a
      % variable
      if contains(string(ib.interfType), "MeshTyingTPFA")
         blocks(b).kind = "variable";
         blocks(b).varName = "pressure";
      else
         blocks(b).kind = "interface";
         blocks(b).varName = ib.coupledVar;
      end
      blocks(b).domainId = [ib.slaveDom, ib.masterDom];
      blocks(b).nDof = ib.nDof;
      blocks(b).solverType = ib.interfType;
      blocks(b).interfId = ib.interfId;
   end

   % Default identity permutation before preconditioner-driven resolution
   pBlk = 1:nTotalBlocks;
   isPermuted = false;
   dofOffsets = [0, cumsum([blocks.nDof])];
   rawDofs = cell(nTotalBlocks, 1);
   for b = 1:nTotalBlocks
      rawDofs{b} = ((dofOffsets(b)+1):dofOffsets(b+1))';
   end
   pDofs = vertcat(rawDofs{:});
   defaultBlocks = blocks;

   % Identify couplings between blocks
   couplings = {};

   % Intra-domain couplings (e.g. Biot, MultiphasePoromechanics, EFEM)
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

   % Inter-domain couplings via interfaces
   for iI = 1:nInterf
      ib = interfBlocks(iI);
      desc = ib.interfType;
      descType = sprintf('interface (%s)', desc);

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
   info.pBlk = pBlk;
   info.pDofs = pDofs;
   info.isPermuted = isPermuted;
   info.blocks = defaultBlocks;
   info.rawBlocks = blocks;
   info.varBlocks = varBlocks;
   info.interfBlocks = interfBlocks;
   info.couplings = couplings;

   % Print summary to console when requested
   if verbose
      fprintf('=== Assembled Matrix Structure (%dx%d blocks) ===\n', nTotalBlocks, nTotalBlocks);
      for b = 1:nTotalBlocks
         blk = defaultBlocks(b);
         if blk.kind == "variable" && ~contains(blk.solverType,"MeshTyingTPFA")
            solvStr = join(blk.solverType, ',');
            fprintf('Block (%d,%d): Domain %d, variable "%s" [DoFs: %d, physics: %s]\n', ...
               b, b, blk.domainId, blk.varName, blk.nDof, solvStr);
         else
            fprintf('Block (%d,%d): Interface (%s) [DoFs: %d, variable: "%s"]\n', ...
               b, b, blk.solverType, blk.nDof, blk.varName);
         end
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

function applyBC(model, bound, t, syst, state)
% Apply Boundary condition to linear system
% BCs is imposed to each physics solver separately
keys = bound.db.keys;
for i = 1 : length(keys)
   cond = bound.getCond(keys{i});
   field = translatePhysic(bound.getPhysics(keys{i}),model);
   type = 'VolumeForce';
   if ~strcmp(cond,'VolumeForce')
      type = bound.getType(keys{i});
   end
   % get id of constrained entities and corresponding BC value
   [bcEnts,bcVals] = getBCVals(keys{i},cond,field,model,bound,t);
   syst.applyBC(type,field,bcEnts,bcVals,state)
end
%
end


function [ents,vals] = getBCVals(bc,cond,field,model,bound,t)
vals = [];
ph = translatePhysic(field,model); % change field tag to query model class
% get entities and values depending on type of BCs
switch cond
   case 'NodeBC'
      ents = bound.getDoF(bc,field);
   case 'ElementBC'
      ents = bound.getDoF(bc,field);
   case 'SurfBC'
      if isFEMBased(model,ph)
         ents = bound.getDoF(bc,field);
         entitiesInfl = bound.getEntitiesInfluence(bc);
         q = bound.getVals(bc, t);
         vals = entitiesInfl*q;
      elseif isFVTPFABased(model,ph)
         ents = bound.getEntities(bc);    % face ID
         vals = bound.getVals(bc, t);
      end
   case 'VolumeForce'
      ents = bound.getDoF(bc,field);
      if isFEMBased(model,ph)
         entitiesInfl = bound.getEntitiesInfluence(bc);
         q = bound.getVals(bc, t);
         vals = entitiesInfl*q;
      elseif isFVTPFABased(model,ph)
         vals = bound.getVals(bc, t);
      end
end
end
function applyBC_test(syst, bcDofs, rhsVal, dirVal, i, bound)
% Apply boundary conditions to each block of the global system
% Penalty method for Dir. cond in case of empty value of rhsVal and dirVal
type = bound.getCond(i);
ph = bound.getPhysics(i);
% get entity type associated to BC
if isFEMBased(bound.model, ph)
    ent = 'node';
elseif isFVTPFABased(bound.model, ph)
    ent = 'elem';
end

% get block ID associated to given bc dofs
phDofs = bound.dof.loc2sub(bcDofs, ent);
switch type
    case 'Dir' % Dirichlet BC
        for i = 1:max(phDofs)
            dof = bcDofs(phDofs == i);
            if isempty(rhsVal) && isempty(dirVal) % penalty method
                nrows = size(syst.blockJ(i,i).block,1);
                maxVal = max(syst.blockJ(i,i).block, [], "all");
                syst.J(nrows*(dof-1) + dof) = maxVal*1.e10;
            else

            end
        end
    case 'Neu'
        for i = 1:max(phDofs)
        end
end






function doflist = getGlobalDofs(list,dofM,dofS,dofIm, flag)
% return dof numbering according to the input flag
switch flag
    case 'master'
        doflist = find(ismember(dofM,list));
    case 'slave'
        doflist = length(dofM) + find(ismember(dofS,list));
    case 'interface'
        doflist = length(dofM) + length(dofS) + find(ismember(dofIm,list));
end


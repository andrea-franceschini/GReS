function doflist = getGlobalDofs(list,dofM,dofS,dofIm,dofIs,flag)
% return dof numbering according to the input flag
switch flag
    case 'master'
        doflist = find(ismember(dofM,list));
    case 'slave'
        doflist = length(dofM) + find(ismember(dofS,list));
    case 'interfaceMaster'
        doflist = length(dofM) + length(dofS) + find(ismember(dofIm,list));
    case 'interfaceSlave'
        doflist = length(dofM) + length(dofS) + length(dofIm) + find(ismember(dofIs,list));
    case 'lagrange'
        doflist = length(dofM) + length(dofS) + length(dofIm) + length(dofIs) + find(ismember(dofIs,list));
end


function [ii,jj,J] = addBlockJ(ii,jj,J,dofI,dofJ,Jblock)
% Add block to sparse matrix definition arrays
[jjAdd,iiAdd] = meshgrid(dofJ,dofI);
% shrink JBlock
% Jblock = Jblock(dofI,dofJ);
ii = [ii; iiAdd(:)];
jj = [jj; jjAdd(:)];
J = [J; Jblock(:)];
end


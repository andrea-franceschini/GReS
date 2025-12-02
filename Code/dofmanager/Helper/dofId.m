function dofOut = dofId(dofIn,nComp)
   dofIn = reshape(dofIn,[],1);        % make sure input is a column vector
   dofOut = nComp*repelem(dofIn,nComp,1)-repmat((nComp-1:-1:0)',size(dofIn,1),1);
end


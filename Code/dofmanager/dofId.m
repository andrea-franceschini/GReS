function dofOut = dofId(dofIn,nComp)
   assert(size(dofIn,2)==1,'Input must be a column vector');
   dofOut = nComp*repelem(dofIn,nComp,1)-repmat((nComp-1:-1:0)',size(dofIn,1),1);
end


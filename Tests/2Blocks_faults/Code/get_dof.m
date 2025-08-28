function dof = get_dof(nodList)
% slaveNodes is a row vector
nodList = reshape(nodList,1,[]);
dof = repelem(3*nodList,3);
dof = (dof + repmat(-2:0,1,numel(nodList)))';
end


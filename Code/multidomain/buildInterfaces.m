function interfaces = buildInterfaces(fileName,domains)

outStruct = readstruct(fileName,AttributeSuffix="");

if isfield(outStruct,"Interface")
  outStruct = outStruct.Interface;
else
  interfaces = [];
  return
end

interfaces = {};

interfNames = fieldnames(outStruct);

k=0;
for i = numel(interfNames)

  % deal with multiple interfaces of the same type
  s = outStruct.(interfNames{i});
  for j = 1:numel(s)

    interf = feval(interfNames{i},...
      k+1,domains,s(j));

    interf.registerInterface(s(j));
    interfaces{end+1} = interf;

    % update interface counter
    k = k+1;
  end
end

end

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
nInt = numel(interfNames);
for i = 1:nInt
  interf = feval(interfNames{i},...
                 domains,outStruct.(interfNames{i}));

  interf.setInterface(outStruct.(interfNames{i}));
  interfaces{end+1} = interf;
end

end

function fID = openReadOnlyFile(fName)
fID = fopen(fName, 'r');
if fID == -1
    error('File %s does not exist in the directory',fName);
end
end


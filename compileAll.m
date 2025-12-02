clc;
clear;
close all;

list = {'Code/read',...
  'Code/write', ...
  'Code/elements',...
  'ThirdPartyLibs/PolyClipper'};

home = pwd;
for folder = list
    fprintf('Compiling MEX files in %s\n', folder{1});
    cd(folder{1});
    compile
    fprintf('\n\n')
    cd(home);
end
clear

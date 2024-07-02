clc;
clear;
close all;

list = {'Code/read', 'Code/write'};

home = pwd;
for folder = list
    fprintf('Compiling MEX files in %s\n', folder{1});
    cd(folder{1});
    compile
    cd(home);
end
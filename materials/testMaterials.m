close all;
clear;
clc;

fileName = 'materials.dat';

mat = Materials(fileName);
tic;
elas = mat.getMaterial('elas');
elas.getStiffnessMatrix()
hypo = mat.getMaterial('hypo');
hypo.getStiffnessMatrix(3)
rock = mat.getMaterial('rock');
rock.getPorosity()
t1 = toc;
fprintf('Time to read %.3f [s]\n', t1);

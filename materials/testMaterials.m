close all;
clear;
clc;

% Setting the input file name
fileName = 'materials.dat';

% Creation of an object of "Materials"
mat = Materials(fileName);
% Starting a timer
tic;

% Calling the "getMaterial" function from the class "Materials" 
elas = mat.getMaterial('elas');

% Calling the function of the class "Elastic" that calculates stiffness
% matrix for elastic materials
D1 = elas.getStiffnessMatrix();

hypopl = mat.getMaterial('hypopl');
% Calling the function of the class "HypoPlastic" that calculates stiffness
% matrix for elastic materials. Z-stress is needed as input data
D2 = hypopl.getStiffnessMatrix(1);

hypoel = mat.getMaterial('hypoel');
% Calling the function of the class "HypoPlastic" that calculates stiffness
% matrix for elastic materials. Z-stress is needed as input data
D3 = hypoel.getStiffnessMatrix(3);

anis = mat.getMaterial('transvel');
D4 = anis.getStiffnessMatrix(3);

% cam = mat.getMaterial('camclay');
% cam.

rock = mat.getMaterial('rock');
% Calling the function of the class "PorousRock" that gets rock
% permeability
rock.getPermeability();

water = mat.getMaterial('fluid');
% Calling the function of the class "Fluid" that gets the fluid weight
water.getWeight();

% Reading the stopwatch timer
t1 = toc;
% Printing t1 = time to read 
fprintf('Time to read %.3f [s]\n', t1);

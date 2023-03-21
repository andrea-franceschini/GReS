close all;
clear all;

% Setting the input file name
fileName = 'materials1.dat';

% Starting a timer
tic;

% Creation and intialization of an object of "Materials"
mat = Materials(fileName);

% Calling the "getMaterial" function from the class "Materials" 
elas = mat.getMaterial(1);

% Calling the function of the class "Elastic" that calculates stiffness
% matrix for elastic materials
D1 = elas.getStiffnessMatrix();

hypoel = mat.getMaterial(2);
% Calling the function of the class "HypoPlastic" that calculates stiffness
% matrix for elastic materials. Z-stress is needed as input data
D3 = hypoel.getStiffnessMatrix(1);
% 
% anis = mat.getMaterial(3);
% D4 = anis.getStiffnessMatrix();
% 
% % cam = mat.getMaterial('camclay');
% % cam.
% 
rock = mat.getMaterial(3);
% Calling the function of the class "PorousRock" that gets rock
% permeability
K = rock.getPermeability();
por = rock.getPorosity();
% 
water = mat.getMaterial(4);
% % Calling the function of the class "Fluid" that gets the fluid weight
gamma = water.getWeight();
beta = water.getCompressibility();

% Reading the stopwatch timer
t1 = toc;
% Printing t1 = time to read 
fprintf('Time to read %.3f [s]\n', t1);

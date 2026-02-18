%% Capillary and relative permeabilities curves generator
% This script generates capillary curves and relative permeability curves
% in a tabular format according to Van Genuchten and Mualem (VGM) models,
% respectively.
%
% CLOSED FORM EQUATION FOR CAPILLARY CURVES (van Genuchten, 1980 and 
% van Genuchten and Nielsen, 1985)
%
% Se = [1 + b^n]^(-m)
%
% where:
% * Se is the effective water saturation (Se=(Sw-Swr)/(1-Swr)), with Sw
%   being the water saturation and Swr the residual water saturation;
% * b = p/p_a with p the fluid pressure and pa the air entry pressure;
% * n and m are two fitting parameters to be found empirically.
%
% m can be expressed in terms of n as follows:
%         m = 1 - 1/n
% This relationship is valid in the range 1.25 < n < 6
%
% ----------------------------------------------------------------------
%
% CLOSED FORM EQUATION FOR RELATIVE PERMEABILITY CURVES (Mualem, 1976)
%
% kr = (1 + b)^(-5/2*m) * [(1 + b)^m - b^m]^2
%
% SIGN CONVENTION: it is assumed that pressure is POSITIVE in the 
% unsaturated zone.
%
close all
clear
%
addpath('./VGParameters');
%
%% Setup the soil database
data = VGMCurves();
%
% Set the van Genuchten/Mualem parameters:
% 1) Use the empirical properties of one soil in the database.
%    The available soils are: Clay, Clayey_loam, Loam, Loamy_sand, Sand,
%    Sandy_clay, Sandy_clayey_loam, Sandy_loam, Silt, Silty_clay,
%    Silty_clayey_loam, Silty_loam.
%    Also provide the water specific weight in the chosen unit.
%
%    Syntax:
%    data.setVGMParameters('Clay','SpecWeight',9.81);
data.setVGMParameters('Sand', 'SpecWeight', 9.81);
% 2) Provide user defined values for n and pa
%
%    Syntax:
%    data.setVGMParameters('User_defined','n', 2, 'pEntry', 10);
% data.setVGMParameters('User_defined','n', 2, 'pEntry', 10);
% data.setVGMParameters('User_defined','n', 3, 'pEntry', 1);
%
%% Compute and print the VGM curves in a tabular format.
% Please specify the following:
% * fileName -> name of the output file for Pc and kr curves, respectively;
% * range -> Exponents of the expected pressure range experienced in the 
%   model (for instance, 10^(-5) < p < 10^3 is [-5, 3]
% * nPoints -> number of points for the approximation of the pc and kr
%   curves.
% fName = ["../OUT/pcCurve_400.dat", "../OUT/krCurve_400.dat"];
fName = ["../OUT/pcCurve_2000.dat", "../OUT/krCurve_2000.dat"];
data.makeCurves('fileName', fName, 'range', [-2,1], 'nPoints', 200);
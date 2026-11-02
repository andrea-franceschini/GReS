close all
clear
xpt = [1 4 10];
ypt = [0.95 0.7 0.1];
dypt = diff(ypt)./diff(xpt);
%
x = [3,-0.1,-1,-6,-8,-2,-10.1,-12];
[y,dy] = interpTable(xpt,ypt,dypt,x);

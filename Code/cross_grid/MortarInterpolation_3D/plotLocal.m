function plotLocal(x,y,z)
% create a surface plot of the RBF interpolated basis function
% within a single slave element
figure(1)
stem3(x, y, z)
grid on
xv = linspace(min(x), max(x), 50);
yv = linspace(min(y), max(y), 50);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);
figure(2)
surf(X, Y, Z);
grid on
shading interp
end

% call: plotLocal(ptsGauss(:,1), ptsGauss(:,2), fSlave)


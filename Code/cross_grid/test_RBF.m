%%
close all; clear

%% IMPORT MESHES
msh1 = Mesh();
msh2 = Mesh();
msh1.importGMSHmesh('meshes/mesh1_quad.msh');
msh2.importGMSHmesh('meshes/mesh1_quad.msh');

%% INTERPOLATE
testFunc = @(x,y,z) sin(2*pi*x).*cos(3*pi*y) + exp(x+y);
fIn = testFunc(msh1.coordinates(:,1), msh1.coordinates(:,2), msh1.coordinates(:,3));
%fOut_test = testFunc(msh2.coordinates(:,1), msh2.coordinates(:,2), msh2.coordinates(:,3));
%fIn = load("sol_terzaghi");
rbf = RBF_interpolation(msh1,msh2,1);
fOut = rbf.interpolate(fIn);




%% PLOT SOLUTION
plotFunction(msh1,'functIn', fIn)
plotFunction(msh2,'functOut', fOut)

% %% COMPUTE ERROR NORM
% % Create an object of the "Elements" class and process the element properties
% elems = Elements(msh2);
% volNod = zeros(msh2.nNodes,1);
% if any(msh2.cellVTKType == 12)
%   N1 = getBasisFinGPoints(elems.hexa);
% end
% for el=1:msh2.nCells
%   top = msh2.cells(el,1:msh2.cellNumVerts(el));
%   if msh2.cellVTKType(el) == 10 % Tetra
%     volNod(top) = volNod(top) + elems.vol(el)/msh2.cellNumVerts(el);
%   elseif msh2.cellVTKType(el) == 12 % Hexa
%     dJWeighed = getDerBasisFAndDet(elems.hexa,el,3);
%     volNod(top) = volNod(top)+ N1'*dJWeighed';
%   end
% end
% 
% %pressure_error
% err2 = (interp.fOut - fOut_test).^2;
% errNorm = sqrt(err2'*volNod);

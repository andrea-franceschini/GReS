function  plotSolution(mesh, foldName, disp, stress)
% Plot stresses and displacements for 2D mesh 
outVTK = VTKOutput(mesh, foldName); % create VTK object

dispData(1,1).name = 'ux';
dispData(2,1).name = 'uy';
dispData(1,1).data = disp(1:2:end);
dispData(2,1).data = disp(2:2:end);

stressData(1,1).name = 'sx';
stressData(2,1).name = 'sy';
stressData(3,1).name = 'tauxy';
stressData(1,1).data = stress(:,1);
stressData(2,1).data = stress(:,2);
stressData(3,1).data = stress(:,3);

outVTK.writeVTKFile(0, [], [], dispData, stressData);
outVTK.finalize();
end


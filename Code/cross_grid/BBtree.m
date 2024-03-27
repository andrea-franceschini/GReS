clear
close all
%% contact search algorithm
msh1 = Mesh();
msh2 = Mesh();
msh1.importGMSHmesh('meshes/mesh1.msh');
msh2.importGMSHmesh('meshes/mesh2.msh');
msh1.computeSurfaceCentroid();
msh2.computeSurfaceCentroid();
cs = ContactSearching(msh1,msh2,8);

%% define arbitrary polygons
n_edges = 4;
l1 = 2*pi*sort(rand(n_edges,1));
l2 = 2*pi*sort(rand(n_edges,1));
xp1 = cos(l1);
yp1 = sin(l1);
xp2 = 1 + cos(l2);
yp2 = 1 + sin(l2);
plot([xp1;xp1(1)],[yp1;yp1(1)])
hold on
plot([xp2;xp2(1)],[yp2;yp2(1)])

%% find bounding boxes
dir = [1 0;
       0 1;
      -1 1;
      1 1];
vals1 = [xp1 yp1]*dir';
[~,imax1] = max(vals1);
[~,imin1] = min(vals1);
vals2 = [xp2 yp2]*dir';
[~,imax2] = max(vals2);
[~,imin2] = min(vals2);
% hold on
% plot(xp(imax),yp(imax), 'ro')
% plot(xp(imin),yp(imin), 'ro')


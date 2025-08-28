clear
close all
clc

% naive algorithm to generate a volume of Pore
% no requested porosity

% input: NX NY NZ for the sphere cloud (tensor grid)
% min and max distance of points in each direction

% after generating the point cloud, every sphere is generated recursively

% the first sphere take a radius such that the radius is smaller than half
% the distance to the other points


% the following sphere follow this criterion: the radius is smaller than:
% half the distance to every other sphere not yet processed 
% OR
% smaller than the distance to other generated spheres (distance + radius)

NX = 50;
NY = 50;
NZ = 50;

Np = 6;

L = 1.5; % size of bounding box

dmin = 0.05;
dmax = 0.3;
x = genPoints1D(dmin,dmax,NX);
y = genPoints1D(dmin,dmax,NY);
z = genPoints1D(dmin,dmax,NZ);

% remap points to 0,1 domain
x = x/max(x);
y = y/max(y);
z = z/max(z);

% tensor grid
[X,Y,Z]= meshgrid(x,y,z);

% discard points in the cloud
id = randi(numel(X),Np,1);

X = X(:);
Y= Y(:);
Z = Z(:);

X = X(id);
Y = Y(id);
Z = Z(id);

scatter3(X(:),Y(:),Z(:))

C = [X, Y, Z];

[C,R] = getRadius(C);

%%
R(1) = 1.05*R(1);
R(3) = 1.1*R(3);

plotSpheres(C,R);


% printfiles
writematrix(R, 'rad.txt', 'Delimiter', ' ');
writematrix(C, 'pts.txt', 'Delimiter', ' ');

%% generate mesh
command = "python Mesh/PoreMech.py";
system(command)


function p = genPoints1D(dmin,dmax,N)
  p = zeros(N,1);
  p(1) = rand*dmin;
  D = dmax - dmin;
  for i=2:N
    p(i) = p(i-1) + dmin + rand*D;
  end

end


function [coord,R] = getRadius(coord)
  % generate spheres recursively
  np = size(coord,1);
  R = zeros(np,1);
  for i = 1:size(coord,1)
    % update coordinate matrix
    cMat = coord([1:i-1 i+1:end],:);
    d = dist(coord(i,:),cMat);
    d = d-R([1:i-1 i+1:end]);
    d(i+1:end) = 0.8*d(i+1:end);
    R(i) = min(d);
  end

  % check that all spheres are connected.
  % if not, remove them!
  isContact = false(np,1);
  for i = 1:size(coord,1)
    cMat = coord([1:i-1 i+1:end],:);
    d = dist(coord(i,:),cMat);
    isContact(i) = any(abs(d-R([1:i-1 i+1:end])-R(i))<1e-3);
  end
  warning('Removed %i floating grains \n',sum(~isContact))
  coord = coord(isContact,:);
  R = R(isContact);  
end

function d = dist(c,coordMat)
  d = sqrt(sum((coordMat-c).^2,2));
end


function plotSpheres(C,R)
% base sphere
[x0,y0,z0] = sphere(15);
figure(1)
hold on
for i = 1:length(R)
  x = C(i,1); y = C(i,2); z = C(i,3);
  surf(R(i)*x0+x,R(i)*y0+y,R(i)*z0+z,'FaceColor',rand(1,3))
end
end


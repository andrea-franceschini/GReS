% check the validity of the contact detection strategy on quadrilateral
% elementsmshM.coordinates
clear 
close all
mshM = Mesh();
mshM.createCartesianGrid(2,1,[-1 1],[-1 1],1,1);
% reordering nodes to anticlockwise order
tmp = mshM.coordinates(4,:);
mshM.coordinates(4,:) = mshM.coordinates(3,:);
mshM.coordinates(3,:) = tmp;
type = 'gauss';
g = Gauss(12,3,2);
elem = Elements(mshM,g);
% deform mesh object changing location of nodes
def = 0.2;
z_m = 0;
M = 5;
def_rig = -0.5+rand;
trans = -def+2*def*rand(4,2);
trans = [trans+def_rig z_m*ones(4,1)];
mshM.coordinates = mshM.coordinates + trans;
ptsRef  = linspace(-1,1,M);
[ys, xs] = meshgrid(ptsRef,ptsRef);
ptsRef = [xs(:), ys(:)];
bf = computeBasisF(elem.quad,ptsRef);
pts = bf*mshM.coordinates;
r = sqrt((max(pts(:,1)) - min(pts(:,1)))^2 + (max(pts(:,2)) - min(pts(:,2)))^2 + (max(pts(:,3)) - min(pts(:,3)))^2);
fiMM = zeros(length(pts),length(pts));
% distance matrix

% build the 4 by 4 matrix
for ii = 1:length(pts)
    dist = (pts(:,1) - pts(ii,1)).^2 + (pts(:,2) - pts(ii,2)).^2 + (pts(:,3) - pts(ii,3)).^2;
    dist = sqrt(dist);
    fiMM(ii,:) = computeRBFentries(dist,type,r);
end

dmat = sqrt((pts(:,1) - pts(:,1)').^2 + (pts(:,2) - pts(:,2)').^2 + (pts(:,3) - pts(:,3)').^2);
fiMMnew = computeRBFentries(dmat,type,r);

% computing weights for support detection
wf = fiMM\bf(:,[1 3]);
w1 = fiMM\ones(length(wf),1);



% wf = fiMM\vals;
% w1 = fiMM\ones(length(pts),1);
% place uniform grid in the reference space of master element
mshS = Mesh();
mshS.createCartesianGrid(2,1,[-10 10],[-10 10],1,1);
% reordering nodes to anticlockwise order
tmp = mshS.coordinates(4,:);
mshS.coordinates(4,:) = mshS.coordinates(3,:);
mshS.coordinates(3,:) = tmp;
% deform mesh object changing location of nodes
def = 0.5;
def_r = 1;
z_s = 0;
def_rig = -def_r+2*def_r*rand;
trans = -def+2*def*rand(4,2);
trans = [trans+def_rig z_s*ones(4,1)];
mshS.coordinates = mshS.coordinates + trans;


N  = 30;
samp = linspace(-1,1,N);
[ys, xs] = meshgrid(samp,samp);
samp = [xs(:), ys(:)];

% compute pts location in the real space of the slave element
bf = computeBasisF(elem.quad,samp);
samp = bf*mshS.coordinates;

fiNM = zeros(length(samp),length(pts));
for ii = 1:length(samp)
    dist = (pts(:,1) - samp(ii,1)).^2 + (pts(:,2) - samp(ii,2)).^2 + (pts(:,3) - samp(ii,3)).^2;
    dist = sqrt(dist);
    fiNM(ii,:) = computeRBFentries(dist,type,r);
end

% evaluating planar functions for support detection
valsOut = (fiNM*wf)./(fiNM*w1);

tol = -1e-6;
% find only positive points
out = find(~all([valsOut>=0, valsOut<=1],2));
%out = find(any(valsOut<tol,2));

%%
% scatter3(pts(:,1),pts(:,2),pts(:,3),"filled")
plot3(mshM.coordinates([1:end 1],1),mshM.coordinates([1:end 1],2),mshM.coordinates([1:end 1],3),'.b-','LineWidth',1,'MarkerSize',30)
hold on
plot3(mshS.coordinates([1:end 1],1),mshS.coordinates([1:end 1],2),mshS.coordinates([1:end 1],3),'.r-','LineWidth',1,'MarkerSize',30)
scatter3(samp(:,1),samp(:,2),samp(:,3),'g','.')
scatter3(samp(:,1),samp(:,2),valsOut,"blue")
scatter3(samp(out,1),samp(out,2),samp(out,3),'.','r')
